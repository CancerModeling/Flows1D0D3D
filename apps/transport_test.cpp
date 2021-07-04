////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <macrocirculation/graph_flow_and_concentration_writer.hpp>
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
#include "macrocirculation/transport.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 1;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  const double t_end = 3.;
  const std::size_t max_iter = 1600000;
  // const std::size_t max_iter = 1;

  const double tau = 1e-4 / 16.;
  const double tau_out = 1e-2;
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  const std::size_t num_macro_edges = 4;
  const std::size_t num_edges_per_segment = 11;

  const mc::PhysicalData physical_data = mc::PhysicalData::set_from_data(400000, 1, 1.028e-2, 9, 1., 1. / num_macro_edges);

  std::cout << mc::calculate_c0(physical_data.G0, physical_data.rho, physical_data.A0) << std::endl;

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  const mc::Point start_point = mc::Point{0, 0, 0};
  const mc::Point end_point = mc::Point{1, 0, 0};
  auto start = graph->create_vertex();
  auto previous_vertex = start;
  for (std::size_t macro_edge_index = 0; macro_edge_index < num_macro_edges; macro_edge_index += 1) {
    auto next_vertex = graph->create_vertex();
    auto vessel = graph->connect(*previous_vertex, *next_vertex, num_edges_per_segment);
    vessel->add_embedding_data(mc::EmbeddingData{
      {mc::convex_combination(start_point, end_point, static_cast<double>(macro_edge_index) / num_macro_edges),
       mc::convex_combination(
         start_point, end_point, static_cast<double>(macro_edge_index + 1) / num_macro_edges)}});
    vessel->add_physical_data(physical_data);
    previous_vertex = next_vertex;
  }
  auto end = previous_vertex;
  end->set_to_free_outflow();

  // set inflow boundary conditions
  // start->set_to_inflow(mc::heart_beat_inflow());
  start->set_to_inflow(mc::heart_beat_inflow(4.));
  //start->set_to_inflow([](double t) { return 1.;});

  graph->finalize_bcs();

  // partition graph
  // TODO: app crashes if not enabled -> fix this!
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // configure solver
  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, false);

  mc::Transport transport(MPI_COMM_WORLD, graph, dof_map_flow, dof_map_transport);

  mc::ExplicitNonlinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  solver.use_ssp_method();

  std::vector<double> u_prev(dof_map_flow->num_dof(), 0);

  interpolate_constant(MPI_COMM_WORLD, *graph, *dof_map_flow, 1., 0, u_prev);
  interpolate_constant(MPI_COMM_WORLD, *graph, *dof_map_flow, 1., 1, u_prev);

  std::vector<mc::Point> points;
  points.reserve(graph->num_edges() * 2);
  std::vector<double> c_vertex_values;
  c_vertex_values.reserve(graph->num_edges() * 2);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  mc::GraphFlowAndConcentrationWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map_flow, dof_map_transport);

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    transport.solve(t, tau, solver.get_solution());
    solver.solve(tau, t);

    t += tau;

    if (it % output_interval == 0) {
      std::cout << "iter " << it << std::endl;
      // if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      //   std::cout << "iter = " << it << ", time = " << solver.get_time() << std::endl;

      // save solution
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport.get_solution(), points, c_vertex_values);

      csv_writer.write(t, solver.get_solution(), transport.get_solution());

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.write(t);
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
