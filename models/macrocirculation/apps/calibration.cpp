////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <argparse.hpp>
#include <chrono>
#include <macrocirculation/communication/mpi.hpp>
#include <memory>
#include <utility>

#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_flow_and_concentration_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/transport.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/windkessel_calibrator.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  argparse::ArgumentParser program("Calibration script for geometries");
  program.add_argument("--input-file")
    .help("path to the input file")
    .default_value(std::string("./data/network-33-vessels.json"));
  program.add_argument("--heart-amplitude")
    .help("the amplitude of a heartbeat")
    .default_value(485.0)
    .action([](const std::string &value) { return std::stod(value); });
  program.add_argument("--tau")
    .help("time step size")
    .default_value(2.5e-4 / 16.)
    .action([](const std::string &value) { return std::stod(value); });
  program.add_argument("--tau-out")
    .help("time step size for the output")
    .default_value(1e-3)
    .action([](const std::string &value) { return std::stod(value); });
  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cout << err.what() << std::endl;
    std::cout << program;
    exit(0);
  }

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append(program.get<std::string>("--input-file"), *graph);

  graph->get_vertex(0)->set_to_inflow(mc::heart_beat_inflow(program.get<double>("--heart-amplitude")));

  // set all vertices to free-outflow
  for (auto v_id : graph->get_vertex_ids()) {
    auto v = graph->get_vertex(v_id);
    if (v->is_windkessel_outflow())
      v->set_to_free_outflow();
  }

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, false);

  const double t_end = 1.;
  const std::size_t max_iter = 160000000;

  const auto tau = program.get<double>("--tau");
  const auto tau_out = program.get<double>("--tau-out");
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);
  std::cout << "tau = " << tau << ", tau_out = " << tau_out << ", output_interval = " << output_interval << std::endl;

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

  mc::GraphFlowAndConcentrationWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map_flow, dof_map_transport);

  mc::WindkesselCalibrator calibrator( graph, true );

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    transport_solver.solve(it * tau, tau, flow_solver.get_solution());
    flow_solver.solve();

    // add total flows
    calibrator.update_flow(flow_solver, tau);

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << flow_solver.get_time() << std::endl;

      csv_writer.write(it * tau, flow_solver.get_solution(), transport_solver.get_solution());

      // double estimate parameters
      calibrator.estimate_parameters();
    }

    // break
    if (flow_solver.get_time() > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
