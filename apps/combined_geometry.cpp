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

#include "macrocirculation/0d_boundary_conditions.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/coupled_explicit_implicit_1d_solver.hpp"
#include "macrocirculation/csv_vessel_tip_writer.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_csv_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/linearized_flow_upwind_evaluator.hpp"
#include "macrocirculation/nonlinear_flow_upwind_evaluator.hpp"
#include "macrocirculation/nonlinear_linear_coupling.hpp"
#include "macrocirculation/petsc/petsc_mat.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/rcr_estimator.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

template<typename SolverType>
void output_tip_values(const mc::GraphStorage &graph, const mc::DofMap &dof_map, const SolverType &solver) {
  for (const auto &v_id : graph.get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    auto &vertex = *graph.get_vertex(v_id);
    auto &edge = *graph.get_edge(vertex.get_edge_neighbors()[0]);
    if (vertex.is_windkessel_outflow() || vertex.is_vessel_tree_outflow()) {
      auto &vertex_dof_indices = dof_map.get_local_dof_map(vertex).dof_indices();

      std::vector<double> vertex_values(vertex_dof_indices.size());
      extract_dof(vertex_dof_indices, solver.get_solution(), vertex_values);

      std::cout << "vertex id = " << vertex.get_id() << " has c = " << vertex_values << std::endl;
    }
  }
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(argv[0], "Combined geometry with explicit implicit solver");
  options.add_options()                                                                                                                                           //
    ("tau", "time step size", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                                                              //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                                                 //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("10"))                                                                  //
    ("calibration", "starts a calibration run", cxxopts::value<bool>()->default_value("false"))                                                                   //
    ("update-interval-transport", "transport is only updated every nth flow step, where n is the update-interval?", cxxopts::value<size_t>()->default_value("5")) //
    ("no-slope-limiter", "Disables the slope limiter", cxxopts::value<bool>()->default_value("false"))                                                            //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  const bool slope_limiter_active = !args["no-slope-limiter"].as<bool>();

  {
    // create_for_node the ascending aorta
    auto graph_nl = std::make_shared<mc::GraphStorage>();

    mc::EmbeddedGraphReader graph_reader;
    graph_reader.append("data/meshes/network-33-vessels-extended.json", *graph_nl);

    auto &v_in = *graph_nl->find_vertex_by_name("cw_in");
    v_in.set_to_inflow_with_fixed_flow(mc::heart_beat_inflow(485.0));

    auto graph_li = std::make_shared<mc::GraphStorage>();
    graph_reader.append("data/meshes/coarse-network-geometry.json", *graph_li);
    graph_reader.set_boundary_data("data/meshes/boundary-combined-geometry-linear-part.json", *graph_li);
    graph_reader.set_boundary_data("data/meshes/boundary-combined-geometry-nonlinear-part.json", *graph_nl);

    mc::set_0d_tree_boundary_conditions(graph_li, "bg_");
    // mc::convert_rcr_to_partitioned_tree_bcs(graph_li);
    // mc::convert_rcr_to_rcl_chain_bcs(graph_li);

    auto coupling = std::make_shared<mc::NonlinearLinearCoupling>(MPI_COMM_WORLD, graph_nl, graph_li);
    coupling->add_coupled_vertices("cw_out_1_1", "bg_141");
    coupling->add_coupled_vertices("cw_out_1_2", "bg_139");
    coupling->add_coupled_vertices("cw_out_1_3", "bg_132");
    coupling->add_coupled_vertices("cw_out_2_1", "bg_135");
    coupling->add_coupled_vertices("cw_out_2_2", "bg_119");

    const bool calibration = args["calibration"].as<bool>();

    if (calibration) {
      for (auto v_id : graph_li->get_vertex_ids()) {
        auto v = graph_li->get_vertex(v_id);
        if (!v->is_leaf())
          continue;
        if (v->is_inflow_with_fixed_flow() || v->is_nonlinear_characteristic_inflow() || v->is_linear_characteristic_inflow())
          continue;
        v->set_to_free_outflow();
        std::cout << "setting  " << v->get_name() << " to free outflow." << std::endl;
      }
      for (auto v_id : graph_nl->get_vertex_ids()) {
        auto v = graph_nl->get_vertex(v_id);
        if (!v->is_leaf())
          continue;
        if (v->is_inflow_with_fixed_flow() || v->is_nonlinear_characteristic_inflow() || v->is_linear_characteristic_inflow())
          continue;
        v->set_to_free_outflow();
        std::cout << "setting  " << v->get_name() << " to free outflow." << std::endl;
      }
    }

    graph_nl->finalize_bcs();
    graph_li->finalize_bcs();

    // mc::naive_mesh_partitioner(*graph_nl, MPI_COMM_WORLD);
    mc::flow_mesh_partitioner(PETSC_COMM_WORLD, *graph_nl, degree);

    // mc::naive_mesh_partitioner(*graph_li, MPI_COMM_WORLD);
    mc::flow_mesh_partitioner(PETSC_COMM_WORLD, *graph_li, degree);

    mc::FlowIntegrator flow_integrator_nl(graph_nl);
    mc::FlowIntegrator flow_integrator_li(graph_li);

    mc::CoupledExplicitImplicit1DSolver solver(MPI_COMM_WORLD, coupling, graph_nl, graph_li, degree, degree);

    const double t_end = args["t-end"].as<double>();
    const std::size_t max_iter = 160000000;

    const auto tau = args["tau"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();

    // const double tau_out = tau;
    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    // configure solver
    solver.setup(tau);

    std::vector<mc::Point> points;
    std::vector<double> p_vertex_values;
    std::vector<double> q_vertex_values;
    std::vector<double> c_vertex_values;

    // initial vessel diameters do not change
    std::vector<double> vessel_A0_li;
    mc::fill_with_vessel_A0(MPI_COMM_WORLD, *graph_li, points, vessel_A0_li);

    // vessels ids do not change, thus we can precalculate them
    std::vector<double> vessel_ids_li;
    mc::fill_with_vessel_id(MPI_COMM_WORLD, *graph_li, points, vessel_ids_li);

    auto dof_map_li = solver.get_implicit_dof_map();
    auto dof_map_nl = solver.get_explicit_dof_map();
    auto solver_li = solver.get_implicit_solver();
    auto solver_nl = solver.get_explicit_solver();
    const auto &u_li = solver_li->get_solution();

    // transport dofmap
    auto dof_map_transport_nl = std::make_shared<mc::DofMap>(*graph_nl);
    auto dof_map_transport_li = std::make_shared<mc::DofMap>(*graph_li);
    mc::DofMap::create_for_transport(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, degree);

    // upwind evaluators:
    auto upwind_evaluator_nl = std::make_shared<mc::NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_nl, dof_map_nl);
    auto variable_upwind_provider_nl = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator_nl, solver.get_explicit_solver());

    // upwind evaluators:
    auto upwind_evaluator_li = std::make_shared<mc::LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_li, dof_map_li);
    auto variable_upwind_provider_li = std::make_shared<mc::UpwindProviderLinearizedFlow>(graph_li, upwind_evaluator_li, solver.get_implicit_solver());

    mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, {variable_upwind_provider_nl, variable_upwind_provider_li}, degree);

    solver_nl->use_ssp_method();

    mc::GraphCSVWriter csv_writer_nl(MPI_COMM_WORLD, "output", "combined_geometry_solution_nl", graph_nl);
    csv_writer_nl.add_setup_data(dof_map_nl, solver_nl->A_component, "a");
    csv_writer_nl.add_setup_data(dof_map_nl, solver_nl->Q_component, "q");
    csv_writer_nl.add_setup_data(dof_map_transport_nl, 0, "c");
    csv_writer_nl.setup();

    mc::GraphCSVWriter csv_writer_li(MPI_COMM_WORLD, "output", "combined_geometry_solution_li", graph_li);
    csv_writer_li.add_setup_data(dof_map_li, solver_li->p_component, "p");
    csv_writer_li.add_setup_data(dof_map_li, solver_li->q_component, "q");
    csv_writer_li.add_setup_data(dof_map_transport_li, 0, "c");
    csv_writer_li.setup();

    mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "combined_geometry_solution");

    mc::CSVVesselTipWriter vessel_tip_writer_li(MPI_COMM_WORLD,
                                                "output", "combined_geometry_solution_tips_li",
                                                graph_li,
                                                {dof_map_li, dof_map_transport_li, transport_solver.get_dof_maps_volume().back()},
                                                {"p", "c", "V"});

    mc::CSVVesselTipWriter vessel_tip_writer_nl(MPI_COMM_WORLD,
                                                "output", "combined_geometry_solution_tips_nl",
                                                graph_nl,
                                                {dof_map_nl, dof_map_transport_nl, transport_solver.get_dof_maps_volume().front()},
                                                {"p", "c", "V"});

    const size_t transport_update_interval = args["update-interval-transport"].as<size_t>();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "updating transport every " << transport_update_interval << " interval" << std::endl;

    const auto begin_t = std::chrono::steady_clock::now();
    double t = 0;
    double t_transport = 0;

    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver.solve(tau, t);
      t += tau;
      if (it % transport_update_interval == 0) {
        upwind_evaluator_nl->init(t_transport + transport_update_interval * tau, solver.get_explicit_solver()->get_solution());
        upwind_evaluator_li->init(t_transport + transport_update_interval * tau, solver.get_implicit_solver()->get_solution());
        transport_solver.solve(transport_update_interval * tau, t_transport + transport_update_interval * tau);
        if (slope_limiter_active)
          transport_solver.apply_slope_limiter(t_transport + transport_update_interval * tau);
        t_transport += transport_update_interval * tau;
        // std::cout << "transport sol-norm" << transport_solver.get_solution().norm2() << std::endl;
        // std::cout << "transport mat-norm" << transport_solver.get_mat().norm1() << std::endl;
        // std::cout << "transport rhs-norm" << transport_solver.get_rhs().norm2() << std::endl;
        // if (transport_solver.get_solution().norm2() == INFINITY)
        // {
        //  std::cout << "noooo!!" << std::endl;
        //  exit(1);
        // }
        // std::cout << "transport sol " << transport_solver.get_solution() << std::endl;
      }

      if (calibration) {
        flow_integrator_li.update_flow(*solver_li, tau);
        flow_integrator_nl.update_flow(*solver_nl, tau);
      }

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", t = " << t << std::endl;

        // std::cout << "  transport sol-norm = " << transport_solver.get_solution().norm2() << std::endl;

        // save solution
        csv_writer_nl.add_data("a", solver_nl->get_solution());
        csv_writer_nl.add_data("q", solver_nl->get_solution());
        csv_writer_nl.add_data("c", transport_solver.get_solution());
        csv_writer_nl.write(t);

        csv_writer_li.add_data("p", solver_li->get_solution());
        csv_writer_li.add_data("q", solver_li->get_solution());
        csv_writer_li.add_data("c", transport_solver.get_solution());
        csv_writer_li.write(t);

        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, solver_li->p_component, u_li, points, p_vertex_values);
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, solver_li->q_component, u_li, points, q_vertex_values);
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_transport_li, 0, transport_solver.get_solution(), points, c_vertex_values);

        pvd_writer.set_points(points);
        pvd_writer.add_vertex_data("p", p_vertex_values);
        pvd_writer.add_vertex_data("q", q_vertex_values);
        pvd_writer.add_vertex_data("vessel_id", vessel_ids_li);
        pvd_writer.add_vertex_data("c", c_vertex_values);
        pvd_writer.add_vertex_data("A", vessel_A0_li);
        pvd_writer.write(t);

        vessel_tip_writer_nl.write(t, {solver_nl->get_solution()}, {transport_solver.get_solution(), transport_solver.get_volumes()});
        vessel_tip_writer_li.write(t, {solver_li->get_solution(), transport_solver.get_solution(), transport_solver.get_volumes()});
      }

      // break
      if (t > t_end + 1e-12)
        break;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

    if (calibration) {
      auto flows_nl = flow_integrator_nl.get_free_outflow_data();
      auto flows_li = flow_integrator_li.get_free_outflow_data();

      std::cout << "nonlinear_flows" << std::endl;
      for (auto &it : flows_nl.flows)
        std::cout << it.first << " " << it.second << std::endl;

      std::cout << "linear_flows" << std::endl;
      for (auto &it : flows_li.flows)
        std::cout << it.first << " " << it.second << std::endl;

      mc::RCREstimator rcr_estimator({graph_nl, graph_li});
      double total_flow = mc::get_total_flow({flows_nl, flows_li});
      auto rcr_parameters_nl = rcr_estimator.estimate_parameters(flows_nl.flows, total_flow);
      auto rcr_parameters_li = rcr_estimator.estimate_parameters(flows_li.flows, total_flow);

      mc::parameters_to_json("data/meshes/boundary-combined-geometry-nonlinear-part.json", rcr_parameters_nl, graph_nl);
      mc::parameters_to_json("data/meshes/boundary-combined-geometry-linear-part.json", rcr_parameters_li, graph_li);
    }
  }

  CHKERRQ(PetscFinalize());
}
