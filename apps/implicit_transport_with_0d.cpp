////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <memory>
#include <utility>

#include "petsc.h"

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/linearized_flow_upwind_evaluator.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/right_hand_side_evaluator.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

void implicit_transport_with_explicit_flow(double tau, double tau_out, double t_end, std::shared_ptr<mc::GraphStorage> graph) {
  const std::size_t max_iter = 1600000;

  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());

  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, 0, true, [](const mc::Vertex& v) -> size_t{
    if (v.is_windkessel_outflow())
      return 1;
    return 0;
  });

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  auto flow_solver = std::make_shared<mc::ExplicitNonlinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver->use_ssp_method();

  auto upwind_evaluator = std::make_shared<mc::NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
  auto variable_upwind_provider = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator, flow_solver);

  mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

  transport_solver.set_inflow_function([](double t) -> double{
    if (t < 1.5)
      return 1.;
    return 0.;
  });

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    flow_solver->solve(tau, t);
    variable_upwind_provider->init(t + tau, flow_solver->get_solution());
    transport_solver.solve(tau, t + tau);

    t += tau;

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", time = " << t << std::endl;

      // save solution
      std::vector<mc::Point> points;
      std::vector<double> c_vertex_values;
      std::vector<double> A_vertex_values;
      std::vector<double> Q_vertex_values;
      std::vector<double> v_vertex_values;
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->A_component, flow_solver->get_solution(), points, A_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->Q_component, flow_solver->get_solution(), points, Q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *variable_upwind_provider, t, points, v_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.add_vertex_data("Q", Q_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("v", v_vertex_values);
      pvd_writer.write(t);

      // output tip values
      for (const auto &v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
        auto &vertex = *graph->get_vertex(v_id);
        auto &edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
        if (vertex.is_windkessel_outflow()) {
          auto &vertex_dof_indices = dof_map_transport->get_local_dof_map(vertex).dof_indices();

          std::vector<double> vertex_values(vertex_dof_indices.size());
          extract_dof(vertex_dof_indices, transport_solver.get_solution(), vertex_values);

          std::cout << "vertex id = " << vertex.get_id() << " has c = " << vertex_values << std::endl;
        }
      }
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
}

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves implicit transport problem"));

  cxxopts::Options options(argv[0], "Implicit transport solver.");
  options.add_options()                                                                                                              //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 8.)))                //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                    //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("6"))                                      //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const double tau = args["tau"].as<double>();
  const double tau_out = args["tau-out"].as<double>();
  const double t_end = args["t-end"].as<double>();

  const std::size_t num_micro_edges = 40;

  // vessel parameters
  //const double vessel_length = 20.5;
  const double vessel_length = 10.;
  const double radius = 0.403;
  const double wall_thickness = 0.067;
  const double elastic_modulus = 400000.0;
  const double density = 1.028e-3;

  auto physical_data_short = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
  auto physical_data_long = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto v0 = graph->create_vertex();
  auto v1 = graph->create_vertex();

  auto edge_0 = graph->connect(*v0, *v1, num_micro_edges);
  edge_0->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1., 0, 0)}});
  edge_0->add_physical_data(physical_data_short);

  v0->set_to_inflow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });
  v1->set_to_windkessel_outflow(1.8, 3870);

  graph->finalize_bcs();

  // partition graph
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  implicit_transport_with_explicit_flow(tau, tau_out, t_end, graph);

  CHKERRQ(PetscFinalize());
}
