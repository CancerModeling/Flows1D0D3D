////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <iomanip>
#include <memory>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_csv_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

namespace test_data {}

std::shared_ptr<mc::GraphStorage> create_3_vessel_network() {
  auto graph = std::make_shared<mc::GraphStorage>();
  auto &v0 = *graph->create_vertex();
  auto &v1 = *graph->create_vertex();
  auto &v2 = *graph->create_vertex();
  auto &v3 = *graph->create_vertex();

  v0.set_name("in");
  v1.set_name("out_1");
  v2.set_name("out_2");
  v3.set_name("out_3");

  const double density = 1.028e-3;

  double vessel_length_0 = 4.0;
  double radius_0 = 1.2;
  double wall_thickness_0 = 0.163;
  double elastic_modulus_0 = 400000.0;
  double gamma_0 = 9;
  size_t number_edges_0 = 10;

  double vessel_length_1 = 2.0;
  double radius_1 = 1.12;
  double wall_thickness_1 = 0.126;
  double elastic_modulus_1 = 400000.0;
  double gamma_1 = 9;
  size_t number_edges_1 = 6;

  double vessel_length_2 = 3.4;
  double radius_2 = 0.62;
  double wall_thickness_2 = 0.08;
  double elastic_modulus_2 = 400000.0;
  double gamma_2 = 9;
  size_t number_edges_2 = 10;

  auto data_0 = mc::PhysicalData::set_from_data(elastic_modulus_0, wall_thickness_0, density, gamma_0, radius_0, vessel_length_0);
  auto data_1 = mc::PhysicalData::set_from_data(elastic_modulus_1, wall_thickness_1, density, gamma_1, radius_1, vessel_length_1);
  auto data_2 = mc::PhysicalData::set_from_data(elastic_modulus_2, wall_thickness_2, density, gamma_2, radius_2, vessel_length_2);

  auto &e0 = *graph->connect(v0, v1, number_edges_0);
  e0.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 0, 0), mc::Point(0, 1, 0)}}));
  auto &e1 = *graph->connect(v1, v2, number_edges_1);
  e1.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 1, 0), mc::Point(1, 1, 0)}}));
  auto &e2 = *graph->connect(v1, v3, number_edges_2);
  e2.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 1, 0), mc::Point(-1, 1, 0)}}));

  e0.add_physical_data(data_0);
  e1.add_physical_data(data_1);
  e2.add_physical_data(data_2);

  v0.set_to_inflow_with_fixed_flow(mc::heart_beat_inflow(485.));

  v2.set_to_windkessel_outflow(1.718414143839568, 0.7369003586207183);
  v3.set_to_windkessel_outflow(6.085065767841017, 0.20809964052427588);

  return graph;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  const double t_end = 1.2;
  const std::size_t max_iter = 160000000;

  const double tau = 5e-5;
  const double tau_out = 1e-1;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  std::cerr << "FIX physical parameters" << std::endl;

  // create_for_node the ascending aorta
  auto graph = create_3_vessel_network();
  graph->finalize_bcs();

  // partition graph
  // TODO: app crashes if not enabled -> fix this!
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // initialize dof map
  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  // configure solver
  mc::ExplicitNonlinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map, degree);
  solver.use_ssp_method();

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;

  mc::GraphCSVWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph);
  csv_writer.add_setup_data(dof_map, solver.A_component, "A");
  csv_writer.add_setup_data(dof_map, solver.Q_component, "Q");
  csv_writer.setup();

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "bifurcation_solution");

  // mc::WindkesselCalibrator calibrator(graph, true);

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {
    solver.solve(tau, t);

    t += tau;

    // calibrator.update_flow(solver, tau);

    if (it % output_interval == 0) {
      std::cout << "iter " << it << std::endl;

      // save solution
      csv_writer.add_data("A", solver.get_solution());
      csv_writer.add_data("Q", solver.get_solution());
      csv_writer.write(t);

      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 0, solver.get_solution(), points, Q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 1, solver.get_solution(), points, A_vertex_values);
      mc::calculate_total_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_total_vertex_values);
      mc::calculate_static_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_static_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("Q", Q_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("p_static", p_static_vertex_values);
      pvd_writer.add_vertex_data("p_total", p_total_vertex_values);
      pvd_writer.write(t);

      // calibrator.estimate_parameters();
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  for (auto e_id : graph->get_edge_ids()) {
    auto &edge = *graph->get_edge(e_id);
    double A, Q;
    solver.evaluate_1d_AQ_values(edge, 0.5, A, Q);
    std::cout << std::scientific << std::setprecision(16) << e_id << " " << A << " " << Q << std::endl;
  }

  const std::vector<double> edge_id_to_A{
    6.5155018925380483e+00,
    6.1024958654010675e+00,
    1.7777985228282513e+00};

  const std::vector<double> edge_id_to_Q{
    4.1595860742808094e+02,
    3.1718725628944367e+02,
    9.1211576554297366e+01};

  for (auto e_id : graph->get_active_edge_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    auto &edge = *graph->get_edge(e_id);
    double A, Q;
    solver.evaluate_1d_AQ_values(edge, 0.5, A, Q);
    std::cout << std::scientific << std::setprecision(16) << e_id << " " << A << " " << Q << std::endl;
  }

  // auto params = calibrator.estimate_parameters();
  // mc::parameters_to_json("params", params, graph);

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;


  MPI_Finalize();
}
