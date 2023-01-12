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
#include <iomanip>

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
#include "macrocirculation/petsc/petsc_mat.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 0;

void implicit_transport_with_implicit_flow(double tau, double tau_out, double t_end, bool apply_slope_limiter, std::shared_ptr<mc::GraphStorage> graph) {
  const std::size_t max_iter = 1600000;

  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  // dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, true);

  auto flow_solver = std::make_shared<mc::ImplicitLinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver->setup(tau);

  auto upwind_evaluator = std::make_shared<mc::LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
  auto variable_upwind_provider = std::make_shared<mc::UpwindProviderLinearizedFlow>(graph, upwind_evaluator, flow_solver);
  // auto variable_upwind_provider = std::make_shared<mc::ConstantUpwindProvider>(100.);

  mc::PetscMat A("lin", *dof_map_flow);
  mc::assemble_matrix_inflow(PETSC_COMM_WORLD, *graph, *dof_map_flow, tau, flow_solver->p_component, flow_solver->q_component, A);
  mc::assemble_matrix_inner_boundaries(PETSC_COMM_WORLD, *graph, *dof_map_flow, tau, flow_solver->p_component, flow_solver->q_component, A);
  mc::assemble_matrix_free_outflow(PETSC_COMM_WORLD, *graph, *dof_map_flow, tau, flow_solver->p_component, flow_solver->q_component, A);
  A.assemble();
  mc::PetscVec rhs(PETSC_COMM_WORLD, "rhs", *dof_map_flow);

  mc::PetscVec u_dst(PETSC_COMM_WORLD, "u_dst", *dof_map_flow);

  mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  std::vector<mc::Point> points;
  std::vector<double> vessel_A0;
  mc::fill_with_vessel_A0(MPI_COMM_WORLD, *graph, points, vessel_A0);

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    flow_solver->solve(tau, t + tau);
    variable_upwind_provider->init(t + tau, flow_solver->get_solution());
    transport_solver.solve(tau, t + tau);
    if (apply_slope_limiter)
      transport_solver.apply_slope_limiter(t + tau);

//    if (t > 0.1)
//    {
//      A.mul(flow_solver->get_solution(), u_dst);
//      rhs.zero();
//      mc::assemble_rhs_inflow(PETSC_COMM_WORLD, *graph, *dof_map_flow, tau, t+tau, flow_solver->p_component, flow_solver->q_component, rhs );
//      u_dst.print();
//      auto& e = *graph->get_edge(0);
//      std::vector< double > p_up (e.num_micro_vertices(), 0);
//      std::vector< double > q_up (e.num_micro_vertices(), 0);
//      upwind_evaluator->get_fluxes_on_macro_edge(t+tau, e, flow_solver->get_solution(), p_up, q_up);
//      double L = mc::linear::get_L(e.get_physical_data());
//      double C = mc::linear::get_C(e.get_physical_data());
//      auto& ldofmap = dof_map_flow->get_local_dof_map(e);
//      std::vector< size_t > dof_indices(ldofmap.num_basis_functions());
//      for (auto meid = 0; meid < e.num_micro_edges(); meid+=1)
//      {
//        double v_q = tau/L * (p_up[meid+1] - p_up[meid]);
//        double v_p = tau/C * (q_up[meid+1] - q_up[meid]);
//
//        std::cout << "edge " << meid << ":" << std::endl;
//        std::cout << std::ios::scientific << std::setprecision(8) << std::endl;
//
//        ldofmap.dof_indices(meid, flow_solver->p_component, dof_indices);
//        double error = std::abs(v_p + rhs.get(dof_indices[0]) - u_dst.get(dof_indices[0]));
//        std::cout << "q " << v_p + rhs.get(dof_indices[0]) << " " << u_dst.get(dof_indices[0]) << " " << error << std::endl;
//
//        if (std::abs(v_p + rhs.get(dof_indices[0]) - u_dst.get(dof_indices[0])) > 1e-10)
//          std::cout << "?" << std::endl;
//        //std::cout << rhs.get(dof_indices[0]) << std::endl;
//
//        ldofmap.dof_indices(meid, flow_solver->q_component, dof_indices);
//        error = std::abs(v_q + rhs.get(dof_indices[0]) - u_dst.get(dof_indices[0]));
//        std::cout << "p " << v_q + rhs.get(dof_indices[0]) << " " << u_dst.get(dof_indices[0]) << " " << error << std::endl;
//        // std::cout << rhs.get(dof_indices[0]) << std::endl;
//
//        if (std::abs(v_q + rhs.get(dof_indices[0]) - u_dst.get(dof_indices[0])) > 1e-10)
//          std::cout << "?" << std::endl;
//      }
//    }

    t += tau;

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", time = " << t << std::endl;

      // save solution
      std::vector<double> c_vertex_values;
      std::vector<double> p_vertex_values;
      std::vector<double> q_vertex_values;
      std::vector<double> v_vertex_values;
      std::vector<double> A_vertex_values;
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->p_component, flow_solver->get_solution(), points, p_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->q_component, flow_solver->get_solution(), points, q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *variable_upwind_provider, t, points, v_vertex_values);

      auto trafo = [](double p, const mc::Edge& e) {
        return e.get_physical_data().A0 + mc::linear::get_C(e.get_physical_data()) * p;
      };
      mc::interpolate_transformation(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->p_component, flow_solver->get_solution(), trafo, points, A_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.add_vertex_data("Q", q_vertex_values);
      pvd_writer.add_vertex_data("p", p_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("v", v_vertex_values);
      pvd_writer.write(t);
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
}

void implicit_transport_with_explicit_flow(double tau, double tau_out, double t_end, bool apply_slope_limiter, std::shared_ptr<mc::GraphStorage> graph) {
  const std::size_t max_iter = 1600000;

  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  auto flow_solver = std::make_shared<mc::ExplicitNonlinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver->use_ssp_method();

  auto upwind_evaluator = std::make_shared<mc::NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
  auto variable_upwind_provider = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator, flow_solver);
  // auto variable_upwind_provider = std::make_shared<mc::ConstantUpwindProvider>(5.);
  // auto variable_upwind_provider = std::make_shared<mc::EmbeddedUpwindProvider>(graph, [](double t, const mc::Point& p){
  //   return /* std::cos( 1 * M_PI * p.x ) */ std::abs(std::cos( M_PI * t / 4. ));
  // });

  mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    flow_solver->solve(tau, t);
    variable_upwind_provider->init(t + tau, flow_solver->get_solution());
    transport_solver.solve(tau, t + tau);
    if (apply_slope_limiter)
      transport_solver.apply_slope_limiter(t + tau);

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
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
}

std::shared_ptr< mc::Edge > add_embedded_edge(mc::GraphStorage& graph,
              mc::Vertex& v1, const mc::Point& p1,
              mc::Vertex& v2, const mc::Point& p2,
              const size_t num_micro_edges,
              bool forward)
{
  std::shared_ptr< mc::Edge > edge;
  if (forward) {
    edge = graph.connect(v1, v2, num_micro_edges);
    edge->add_embedding_data({{p1, p2}});
  } else {
    edge = graph.connect(v2, v1, num_micro_edges);
    edge->add_embedding_data({{p2, p1}});
  }
  return edge;
}

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves implicit transport problem"));

  cxxopts::Options options(argv[0], "Implicit transport solver.");
  options.add_options()                                                                                                              //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                    //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("6"))                                      //
    ("no-upper-vessel", "Disables the upper vessel at the bifurcation", cxxopts::value<bool>()->default_value("false"))              //
    ("no-lower-vessel", "Disables the lower vessel at the bifurcation", cxxopts::value<bool>()->default_value("false"))              //
    ("flip-left", "", cxxopts::value<bool>()->default_value("false"))              //
    ("flip-right", "", cxxopts::value<bool>()->default_value("false"))              //
    ("flip-upper", "", cxxopts::value<bool>()->default_value("false"))              //
    ("flip-lower", "", cxxopts::value<bool>()->default_value("false"))              //
    ("apply-slope-limiter", "Applies a slope limiter the the transport", cxxopts::value<bool>()->default_value("false"))              //
    ("explicit-flow", "Enables an explicit flow solver instead of the implicit one", cxxopts::value<bool>()->default_value("false")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const double tau = args["tau"].as<double>();
  std::cout << "tau = " << tau << std::endl;
  const double tau_out = args["tau-out"].as<double>();
  const double t_end = args["t-end"].as<double>();

  const bool no_lower_vessel = args["no-lower-vessel"].as<bool>();
  const bool no_upper_vessel = args["no-upper-vessel"].as<bool>();

  const std::size_t num_micro_edges = 20;

  // vessel parameters
  //const double vessel_length = 20.5;
  const double vessel_length = 5.25;
  const double radius = 0.403;
  const double wall_thickness = 0.067;
  const double elastic_modulus = 400000.0;
  const double density = 1.028e-3;

  auto physical_data_short = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
  auto physical_data_long = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);

  const bool flip_left = args["flip-left"].as< bool >();
  const bool flip_right = args["flip-right"].as< bool >();
  const bool flip_upper = args["flip-upper"].as< bool >();
  const bool flip_lower = args["flip-lower"].as< bool >();

  const bool apply_slope_limiter = args["apply-slope-limiter"].as< bool >();

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto v0 = graph->create_vertex();
  auto v0m1 = graph->create_vertex();
  auto v1 = graph->create_vertex();
  auto v2 = graph->create_vertex();

  auto edge_0 = graph->connect(*v0, *v0m1, num_micro_edges);
  edge_0->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(0.25, 0, 0)}});
  edge_0->add_physical_data(physical_data_short);

  auto edge_1 = add_embedded_edge(*graph, *v0m1, mc::Point(0.25, 0, 0), *v1, mc::Point(0.5, 0, 0), num_micro_edges, !flip_left);
  if (flip_left)
    std::cout << "flips left edge" << std::endl;
  edge_1->add_physical_data(physical_data_short);

  auto edge_2 = add_embedded_edge(*graph, *v1,mc::Point(0.5, 0, 0), *v2, mc::Point(1., 0, 0), num_micro_edges, !flip_right);
  if (flip_right)
    std::cout << "flips right edge" << std::endl;
  edge_2->add_physical_data(physical_data_long);

  if (!no_upper_vessel) {
    auto v3 = graph->create_vertex();
    auto edge_3 = add_embedded_edge(*graph, *v1,mc::Point(0.5, 0, 0), *v3, mc::Point(0.5, 0.5, 0), num_micro_edges, !flip_upper);
    if (flip_upper)
      std::cout << "flips upper edge" << std::endl;
    edge_3->add_physical_data(physical_data_long);
    v3->set_to_free_outflow();
  }

  if (!no_lower_vessel) {
    auto v4 = graph->create_vertex();
    auto edge_4 = add_embedded_edge(*graph, *v1, mc::Point(0.5, 0, 0), *v4, mc::Point(0.5, -0.5, 0), num_micro_edges, !flip_lower);
    if (flip_lower)
      std::cout << "flips lower edge" << std::endl;
    edge_4->add_physical_data(physical_data_long);
    v4->set_to_free_outflow();
  }

  v0->set_to_inflow_with_fixed_flow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });
  v2->set_to_free_outflow();

  graph->finalize_bcs();

  // partition graph
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  const bool explicit_flow = args["explicit-flow"].as<bool>();

  if (explicit_flow) {
    implicit_transport_with_explicit_flow(tau, tau_out, t_end, apply_slope_limiter, graph);
  } else {
    implicit_transport_with_implicit_flow(tau, tau_out, t_end, apply_slope_limiter, graph);
  }

  CHKERRQ(PetscFinalize());
}
