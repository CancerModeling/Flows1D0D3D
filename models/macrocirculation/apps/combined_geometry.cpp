////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <macrocirculation/communication/mpi.hpp>
#include <macrocirculation/dof_map.hpp>
#include <macrocirculation/explicit_nonlinear_flow_solver.hpp>
#include <macrocirculation/graph_csv_writer.hpp>
#include <macrocirculation/graph_flow_and_concentration_writer.hpp>
#include <macrocirculation/graph_pvd_writer.hpp>
#include <macrocirculation/interpolate_to_vertices.hpp>
#include <macrocirculation/quantities_of_interest.hpp>
#include <macrocirculation/transport.hpp>
#include <memory>

#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

void connect(mc::GraphStorage &graph, mc::Vertex &start, mc::Vertex &end, const size_t num_micro_edges, const std::vector<mc::PhysicalData> &physical_data) {
  const auto num_macro_edges = physical_data.size();

  mc::Vertex *previous_vertex = &start;
  for (std::size_t macro_edge_index = 0; macro_edge_index < num_macro_edges - 1; macro_edge_index += 1) {
    auto next_vertex = graph.create_vertex();
    auto vessel = graph.connect(*previous_vertex, *next_vertex, num_micro_edges);
    vessel->add_physical_data(physical_data[macro_edge_index]);
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "added vessel with id = " << vessel->get_id() << "(length = " << physical_data[macro_edge_index].length << ", #edges = " << num_micro_edges << ")" << std::endl;
    previous_vertex = next_vertex.get();
  }

  // connect to end point
  auto vessel = graph.connect(*previous_vertex, end, num_micro_edges);
  vessel->add_physical_data(physical_data.back());

  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "added vessel with id = " << vessel->get_id() << "(length = " << physical_data.back().length << ", #edges = " << num_micro_edges << ")" << std::endl;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  // graph_reader.append("data/coarse-network-geometry.json", *graph);
  graph_reader.append("data/network-33-vessels-with-connection.json", *graph);

  auto &v_in = *graph->get_vertex(0);
  auto &v1 = *graph->get_vertex(30);
  auto &v2 = *graph->get_vertex(31);

  v_in.set_to_inflow(mc::heart_beat_inflow(485.0));

  graph_reader.append("data/coarse-network-geometry.json", *graph);

  auto create_physical_data = [](double radius, double vessel_length) -> mc::PhysicalData {
    const double wall_thickness = 0.1;
    const double elastic_modulus = 1.3e6 / 100;
    const double A0 = std::pow(radius, 2) * M_PI;
    const double G0 = 4.0 / 3.0 * std::sqrt(M_PI) * elastic_modulus * wall_thickness / std::sqrt(A0);
    const double rho = 1.028e-3;
    const double viscosity = 4.24e-2; // 4.5 mPa / s
    const double gamma = 9;
    return {elastic_modulus, G0, A0, rho, vessel_length, viscosity, gamma, radius};
  };

  auto interpolate_physical_data = [create_physical_data](size_t num_macro_edges, double start_radius, double end_radius, double vessel_length) {
    std::vector<mc::PhysicalData> physical_data;
    for (int i = 1; i <= num_macro_edges; i += 1) {
      auto r = start_radius * double(num_macro_edges + 1 - i) / double(num_macro_edges + 1) + end_radius * double(i) / double(num_macro_edges + 1);
      physical_data.push_back(create_physical_data(r, vessel_length / num_macro_edges));
    }
    return physical_data;
  };

  // connect both geometries
  const double vessel_length = 10;
  const size_t num_micro_edges = 60;
  const size_t num_macro_edges = 10;

  {
    auto inflow_vertices = graph->find_embedded_vertices({9.093333333333334, 9.173333333333334, 8.053333333333335});
    // std::vector< mc::Vertex* > inflow_vertices = { graph->create_vertex().get() };
    // inflow_vertices[0]->set_to_free_outflow();
    if (inflow_vertices.size() != 1) {
      std::cerr << "expected to find a single vertex, not " << std::to_string(inflow_vertices.size()) << " vertices" << std::endl;
      exit(-1);
    }

    auto physical_data = interpolate_physical_data(num_macro_edges, 0.24, 0.04, vessel_length);
    connect(*graph, v1, *inflow_vertices[0], std::ceil(num_micro_edges / num_macro_edges), physical_data);

    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "finished inter-network connection" << std::endl;
  }

  {
    auto inflow_vertices = graph->find_embedded_vertices({2.9333333333333336, 9.973333333333334, 10.933333333333334});
    // std::vector< mc::Vertex* > inflow_vertices = { graph->create_vertex().get() };
    // inflow_vertices[0]->set_to_free_outflow();
    if (inflow_vertices.size() != 1) {
      std::cerr << "expected to find a single vertex, not " << std::to_string(inflow_vertices.size()) << " vertices" << std::endl;
      exit(-1);
    }

    auto physical_data = interpolate_physical_data(num_macro_edges, 0.24, 0.026, vessel_length);
    connect(*graph, v1, *inflow_vertices[0], std::ceil(num_micro_edges / num_macro_edges), physical_data);

    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "finished inter-network connection" << std::endl;
  }

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, false);

  const double t_end = 100;
  const std::size_t max_iter = 1600000000;

  const double tau = 2.5e-4 / 32 / 4 / 2;
  const double tau_out = 1e-2;
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> flow_solver(MPI_COMM_WORLD, graph, dof_map_flow);
  flow_solver.set_tau(tau);
  flow_solver.use_ssp_method();

  mc::Transport transport_solver(MPI_COMM_WORLD, graph, dof_map_flow, dof_map_transport);

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;
  std::vector<double> c_vertex_values;
  std::vector<double> vessel_ids;

  // vessels ids do not change, thus we can precalculate them
  mc::fill_with_vessel_id(MPI_COMM_WORLD, *graph, points, vessel_ids);

  mc::GraphFlowAndConcentrationWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map_flow, dof_map_transport);
  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "combined_geometry_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    transport_solver.solve(it * tau, tau, flow_solver.get_solution());
    flow_solver.solve();

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << flow_solver.get_time() << std::endl;

      // save solution
      csv_writer.write(flow_solver.get_time(), flow_solver.get_solution(), transport_solver.get_solution());

      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, 0, flow_solver.get_solution(), points, Q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, 1, flow_solver.get_solution(), points, A_vertex_values);
      mc::calculate_total_pressure(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver.get_solution(), points, p_total_vertex_values);
      mc::calculate_static_pressure(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver.get_solution(), points, p_static_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("Q", Q_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("p_static", p_static_vertex_values);
      pvd_writer.add_vertex_data("p_total", p_total_vertex_values);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.add_vertex_data("vessel_id", vessel_ids);
      pvd_writer.write(flow_solver.get_time());
    }

    // break
    if (flow_solver.get_time() > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
