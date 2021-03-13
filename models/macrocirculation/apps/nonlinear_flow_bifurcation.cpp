////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <memory>

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

  const double t_end = 4e-1;
  const std::size_t max_iter = 160000000;

  const double tau = 2.5e-4 / 32 / 4;
  const double tau_out = 1e-3;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  const std::size_t num_edges_per_segment = 80;

  // we create data for celiac ii
  mc::PhysicalData main_vessel_data{
    1706.7e2, // 1706.7 hPa
    0.13,     // 0.13 cm^2
    1.028,    // 1.028 kg/cm^3,
    1.        // length
  };

  mc::PhysicalData other_vessel_data{
    1706.7e2 / std::sqrt(2), // [hPa]
    0.13 / 2,                // [cm^2]
    1.028,                   // [kg/cm^3]
    std::sqrt(2) / 2.        // [cm]
  };

  // create the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  // create vertices:
  auto start = graph->create_vertex();
  auto midpoint1 = graph->create_vertex();
  auto upper_point = graph->create_vertex();
  auto lower_point = graph->create_vertex();
  auto midpoint2 = graph->create_vertex();
  auto endpoint = graph->create_vertex();
  // connect vertices:
  auto e1 = graph->connect(*start, *midpoint1);
  e1->add_physical_data(main_vessel_data);
  auto e2_1 = graph->connect(*midpoint1, *upper_point);
  e2_1->add_physical_data(other_vessel_data);
  auto e2_2 = graph->connect(*midpoint1, *lower_point);
  e2_2->add_physical_data(other_vessel_data);
  auto e3_1 = graph->connect(*lower_point, *midpoint2);
  e3_1->add_physical_data(other_vessel_data);
  auto e3_2 = graph->connect(*upper_point, *midpoint2);
  e3_2->add_physical_data(other_vessel_data);
  auto e4 = graph->connect(*midpoint2, *endpoint);
  e4->add_physical_data(main_vessel_data);
  // to embed the graph in 3D, we define and assign coordinates to the edges:
  auto start_coord = mc::Point(0, 0, 0);
  auto midpoint1_coord = mc::Point(1, 0, 0);
  auto upper_point_coord = mc::Point(2, 1, 0);
  auto lower_point_coord = mc::Point(2, -1, 0);
  auto midpoint2_coord = mc::Point(3, 0, 0);
  auto endpoint_coord = mc::Point(4, 0, 0);
  e1->add_embedding_data({ {start_coord, midpoint1_coord } });
  e2_1->add_embedding_data({ {midpoint1_coord, upper_point_coord } });
  e2_2->add_embedding_data({ {midpoint1_coord, lower_point_coord } });
  e3_1->add_embedding_data({ {upper_point_coord, midpoint2_coord } });
  e3_2->add_embedding_data({ {lower_point_coord, midpoint2_coord } });
  e4->add_embedding_data({ {midpoint2_coord ,endpoint_coord } });

  // set inflow boundary conditions
  start->set_to_inflow(mc::heart_beat_inflow());

  // partition graph
  // TODO: app crashes if not enabled -> fix this!
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // initialize dof map
  auto dof_map = std::make_shared<mc::DofMap>(graph->num_edges());
  fill(MPI_COMM_WORLD, *graph, *dof_map, 2, degree, num_edges_per_segment);

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
  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "bifurcation_solution");

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
