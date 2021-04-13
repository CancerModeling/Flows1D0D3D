////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <macrocirculation/dof_map.hpp>
#include <macrocirculation/explicit_nonlinear_flow_solver.hpp>
#include <macrocirculation/graph_csv_writer.hpp>
#include <macrocirculation/graph_pvd_writer.hpp>
#include <macrocirculation/interpolate_to_vertices.hpp>
#include <macrocirculation/quantities_of_interest.hpp>
#include <memory>

#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append("coarse-network-geometry.json", *graph);

  // auto inflow_vertices = graph->find_embedded_vertices({ 9.093333333333334, 9.173333333333334, 8.053333333333335 });
  std::vector<mc::Point> inflow_points = {
    {9.093333333333334, 9.173333333333334, 8.053333333333335},
    {2.9333333333333336, 9.973333333333334, 10.933333333333334}
  };
  for (auto &p : inflow_points) {
    auto inflow_vertices = graph->find_embedded_vertices(p );
    if (inflow_vertices.size() != 1) {
      std::cerr << "expected to find a single vertex, not " << std::to_string(inflow_vertices.size()) << " vertices" << std::endl;
      exit(-1);
    }
    inflow_vertices[0]->set_to_inflow(mc::heart_beat_inflow(4.85 / 8.));
  }

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map = std::make_shared<mc::DofMap>(graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  /*
  std::vector< mc::Point > points;
  std::vector< double > u_vertex_values;

  std::vector< double > u(dof_map->num_dof(), 0);
  mc::set_to_A0(MPI_COMM_WORLD, *graph, *dof_map, u);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "breast_geometry_solution");
  mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 1, u, points, u_vertex_values);

  pvd_writer.set_points(points);
  pvd_writer.add_vertex_data("A", u_vertex_values);
  pvd_writer.write(0);
   */

  const double t_end = 2.;
  const std::size_t max_iter = 160000000;

  const double tau = 2.5e-4 / 32 / 4 / 2;
  const double tau_out = 1e-3;
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> solver(MPI_COMM_WORLD, graph, dof_map);
  solver.set_tau(tau);
  solver.use_ssp_method();

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;

  mc::GraphCSVWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map, {"Q", "A"});
  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "breast_geometry_solution");

  {
    mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 0, solver.get_solution(), points, Q_vertex_values);
    mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 1, solver.get_solution(), points, A_vertex_values);
    mc::calculate_total_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_total_vertex_values);
    mc::calculate_static_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_static_vertex_values);

    pvd_writer.set_points(points);
    pvd_writer.add_vertex_data("Q", Q_vertex_values);
    pvd_writer.add_vertex_data("A", A_vertex_values);
    pvd_writer.add_vertex_data("p_static", p_static_vertex_values);
    pvd_writer.add_vertex_data("p_total", p_total_vertex_values);
    pvd_writer.write(solver.get_time());
  }

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    solver.solve();

    if (it % output_interval == 0) {
      std::cout << "iter " << it << std::endl;

      // save solution
      csv_writer.write(solver.get_time(), solver.get_solution());

      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 0, solver.get_solution(), points, Q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 1, solver.get_solution(), points, A_vertex_values);
      mc::calculate_total_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_total_vertex_values);
      mc::calculate_static_pressure(MPI_COMM_WORLD, *graph, *dof_map, solver.get_solution(), points, p_static_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("Q", Q_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("p_static", p_static_vertex_values);
      pvd_writer.add_vertex_data("p_total", p_total_vertex_values);
      pvd_writer.write(solver.get_time());
    }

    // break
    if (solver.get_time() > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
