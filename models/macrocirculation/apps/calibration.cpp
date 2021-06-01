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

  cxxopts::Options options(argv[0], "Abstract 33 vessel geometry");
  options.add_options()                                                                                                      //
    ("input-file", "path to the input file", cxxopts::value<std::string>()->default_value("./data/network-33-vessels.json")) //
    ("heart-amplitude", "the amplitude of a heartbeat", cxxopts::value<double>()->default_value("485.0"))                    //
    ("tau", "time step size", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                         //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-3"))                            //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if (!args.unmatched().empty())
  {
    std::cout << "The following arguments were unmatched: " << std::endl;
    for (auto &it : args.unmatched())
      std::cout << " " << it;
    std::cout << "\nAre they part of petsc or a different auxillary library?" << std::endl;
  }

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append(args["input-file"].as<std::string>(), *graph);

  graph->get_vertex(0)->set_to_inflow(mc::heart_beat_inflow(args["heart-amplitude"].as<double>()));

  // set all vertices to free-outflow
  for (auto v_id : graph->get_vertex_ids()) {
    auto v = graph->get_vertex(v_id);
    if (v->is_windkessel_outflow())
      v->set_to_free_outflow();
  }

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  const double t_end = 1.;
  const std::size_t max_iter = 160000000;

  const auto tau = args["tau"].as<double>();
  const auto tau_out = args["tau-out"].as<double>();
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);
  std::cout << "tau = " << tau << ", tau_out = " << tau_out << ", output_interval = " << output_interval << std::endl;

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> flow_solver(MPI_COMM_WORLD, graph, dof_map_flow);
  flow_solver.set_tau(tau);
  flow_solver.use_ssp_method();

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;

  mc::WindkesselCalibrator calibrator(graph, true);

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

  mc::parameters_to_json("data/boundary-data-33-vessels.json", calibrator.estimate_parameters(), graph);

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
