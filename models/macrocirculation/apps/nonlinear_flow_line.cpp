////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
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

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  const double t_end = 0.4;
  const std::size_t max_iter = 1600000;
  // const std::size_t max_iter = 1;

  const double tau = 2.5e-4 / 4 / 2 / 6;
  const double tau_out = 1e-3;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  const std::size_t num_macro_edges = 4 * 6;
  const std::size_t num_edges_per_segment = 44;

  const mc::PhysicalData physical_data = {
    592.4e2,               // 592.4 10^2 Pa,
    6.97,                  // 6.97 cm^2,
    1.028,                 // 1.028 kg/cm^3,
    1.0 / num_macro_edges // 4 cm
  };

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  const mc::Point start_point = mc::Point{0, 0, 0};
  const mc::Point end_point = mc::Point{4, 0, 0};
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
  start->set_to_inflow(mc::heart_beat_inflow());

  // partition graph
  // TODO: app crashes if not enabled -> fix this!
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // configure solver
  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  mc::ExplicitNonlinearFlowSolver<degree> solver(MPI_COMM_WORLD, graph, dof_map);
  solver.set_tau(tau);
  solver.use_ssp_method();

  std::vector<mc::Point> points;
  points.reserve(graph->num_edges() * 2);
  std::vector<double> Q_vertex_values;
  Q_vertex_values.reserve(graph->num_edges() * 2);
  std::vector<double> A_vertex_values;
  A_vertex_values.reserve(graph->num_edges() * 2);
  std::vector<double> p_total_vertex_values;
  p_total_vertex_values.reserve(graph->num_edges() * 2);
  std::vector<double> p_static_vertex_values;
  p_static_vertex_values.reserve(graph->num_edges() * 2);

  mc::GraphCSVWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map, {"Q", "A"});
  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "line_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    solver.solve();

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", time = " << solver.get_time() << std::endl;

      csv_writer.write(solver.get_time(), solver.get_solution());

      // save solution
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
