////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <macrocirculation/communication/mpi.hpp>
#include <memory>
#include <utility>

#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/windkessel_calibrator.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  cxxopts::Options options(argv[0], "Abstract 33 vessel geometry");
  options.add_options()                                                                                                              //
    ("input-file", "path to the input file", cxxopts::value<std::string>()->default_value("./data/meshes/network-33-vessels.json"))         //
    ("output-file", "path to the output file", cxxopts::value<std::string>()->default_value("./data/meshes/boundary-data-33-vessels.json")) //
    ("inflow-vertex-name", "the name of the inflow vertex", cxxopts::value<std::string>()->default_value("cw_in"))                   //
    ("heart-amplitude", "the amplitude of a heartbeat", cxxopts::value<double>()->default_value("485.0"))                            //
    ("heart-period", "the period of one heartbeat", cxxopts::value<double>()->default_value("1."))                                   //
    ("heart-systole-period", "the period of one heartbeat", cxxopts::value<double>()->default_value("0.3"))                          //
    ("tau", "time step size", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                                 //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                    //
    ("t-end", "Time when our simulation ends", cxxopts::value<double>()->default_value("1."))                                        //
    ("verbose", "verbose output", cxxopts::value<bool>()->default_value("false"))                                                    //
    ("h,help", "print usage");
  // options.allow_unrecognised_options(); // for petsc, but we do not use petsc here :P
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append(args["input-file"].as<std::string>(), *graph);

  auto heart_amplitude = args["heart-amplitude"].as<double>();
  auto heart_period = args["heart-period"].as<double>();
  auto heart_systole_period = args["heart-systole-period"].as<double>();
  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "heart amplitude = " << heart_amplitude << " period = " << heart_period << " systole-period = " << heart_systole_period << std::endl;
  auto heart = mc::heart_beat_inflow(heart_amplitude, heart_period, heart_systole_period);
  graph->find_vertex_by_name(args["inflow-vertex-name"].as<std::string>())->set_to_inflow(heart);

  // set all vertices to free-outflow
  for (auto v_id : graph->get_vertex_ids()) {
    auto v = graph->get_vertex(v_id);
    if (v->is_windkessel_outflow())
      v->set_to_free_outflow();
  }

  graph->finalize_bcs();

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  const double t_end = args["t-end"].as<double>();
  const std::size_t max_iter = 160000000;

  const auto tau = args["tau"].as<double>();
  const auto tau_out = args["tau-out"].as<double>();
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);
  std::cout << "tau = " << tau << ", tau_out = " << tau_out << ", output_interval = " << output_interval << std::endl;

  // configure solver
  mc::ExplicitNonlinearFlowSolver flow_solver(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver.set_tau(tau);
  flow_solver.use_ssp_method();

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;

  mc::WindkesselCalibrator calibrator(graph, args["verbose"].as<bool>());

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    flow_solver.solve();

    // add total flows
    calibrator.update_flow(flow_solver, tau);

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << flow_solver.get_time() << std::endl;

      // double estimate parameters
      calibrator.estimate_parameters_local();
    }

    // break
    if (flow_solver.get_time() > t_end + 1e-12)
      break;
  }

  mc::parameters_to_json(args["output-file"].as<std::string>(), calibrator.estimate_parameters(), graph);

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
