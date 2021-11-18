////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "linearized_heart_to_breast_1d_solver.hpp"

#include "0d_boundary_conditions.hpp"
#include "communication/mpi.hpp"
#include "csv_vessel_tip_writer.hpp"
#include "dof_map.hpp"
#include "embedded_graph_reader.hpp"
#include "graph_csv_writer.hpp"
#include "graph_partitioner.hpp"
#include "graph_pvd_writer.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "implicit_transport_solver.hpp"
#include "interpolate_to_vertices.hpp"
#include "linearized_flow_upwind_evaluator.hpp"
#include "vessel_formulas.hpp"

// TODO:  this has to be included because of the get_vessel_tip_dof_values function -> refactor it by moving it somewhere else.
#include "heart_to_breast_1d_solver.hpp"

namespace macrocirculation {

LinearizedHeartToBreast1DSolver::LinearizedHeartToBreast1DSolver(MPI_Comm comm)
    : d_comm(comm),
      d_degree(2),
      d_is_setup(false),
      graph{std::make_shared<GraphStorage>()},
      d_tau_flow{NAN},
      solver{nullptr} {}

void LinearizedHeartToBreast1DSolver::set_path_inflow_pressures(const std::string &path) {
  if (d_is_setup)
    throw std::runtime_error("HeartToBreast1DSolver::set_path_inflow_pressures: cannot be called after setup.");

  path_inflow_pressures = path;
}

void LinearizedHeartToBreast1DSolver::set_path_geometry(const std::string &path) {
  if (d_is_setup)
    throw std::runtime_error("HeartToBreast1DSolver::set_path_geometry: cannot be called after setup.");

  path_linear_geometry = path;
}

void LinearizedHeartToBreast1DSolver::setup(size_t degree, double tau) {
  d_is_setup = true;
  setup_graphs();
  setup_solver(degree, tau);
  setup_output();
}

void LinearizedHeartToBreast1DSolver::solve_flow(double tau, double t) {
  if (std::abs(tau - d_tau_flow) > 1e-12)
    throw std::runtime_error("changing time step width not supported yet.");

  solver->solve(tau, t);
}

double LinearizedHeartToBreast1DSolver::solve_flow(double tau, double t, size_t num_iter) {
  for (size_t idx = 0; idx < num_iter; idx += 1)
  {
    solve_flow(tau, t);
    t += tau;
  }
  return t;
}

void LinearizedHeartToBreast1DSolver::solve_transport(double tau, double t) {
  upwind_evaluator->init(t, solver->get_solution());
  transport_solver->solve(tau, t);
}

void LinearizedHeartToBreast1DSolver::apply_slope_limiter_transport(double t) {
  transport_solver->apply_slope_limiter(t);
}

void LinearizedHeartToBreast1DSolver::write_output(double t) {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto &dof_map = solver->get_dof_map();

  auto &dof_map_transport = *transport_solver->get_dof_maps_transport().at(0);

  csv_writer->add_data("p", solver->get_solution());
  csv_writer->add_data("q", solver->get_solution());
  csv_writer->add_data("c", transport_solver->get_solution());
  csv_writer->write(t);

  interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map, solver->p_component, solver->get_solution(), points, p_vertex_values);
  interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map, solver->q_component, solver->get_solution(), points, q_vertex_values);
  interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map_transport, 0, transport_solver->get_solution(), points, c_vertex_values);

  graph_pvd_writer->set_points(points);
  graph_pvd_writer->add_vertex_data("p", p_vertex_values);
  graph_pvd_writer->add_vertex_data("q", q_vertex_values);
  graph_pvd_writer->add_vertex_data("c", c_vertex_values);
  graph_pvd_writer->add_vertex_data("vessel_id", vessel_ids);
  graph_pvd_writer->add_vertex_data("r", vessel_radii);
  graph_pvd_writer->add_vertex_data("A", vessel_A);
  graph_pvd_writer->write(t);

  vessel_tip_writer->write(t, {solver->get_solution(), transport_solver->get_solution(), transport_solver->get_volumes()});
}

void LinearizedHeartToBreast1DSolver::set_output_folder(std::string output_dir) {
  output_folder_name = std::move(output_dir);
}

std::vector<VesselTipCurrentCouplingData> LinearizedHeartToBreast1DSolver::get_vessel_tip_pressures() {
  auto &u_flow = solver->get_solution();
  auto pressure_values = get_vessel_tip_dof_values(d_comm, *graph, *d_dof_map, u_flow);

  auto &dof_map_transport = *transport_solver->get_dof_maps_transport().back();
  auto &u_transport = transport_solver->get_solution();
  auto concentration_values = get_vessel_tip_dof_values(d_comm, *graph, dof_map_transport, u_transport);

  std::vector<VesselTipCurrentCouplingData> results;

  for (auto v_id : graph->get_vertex_ids()) {
    auto &v = *graph->get_vertex(v_id);

    if (v.is_vessel_tree_outflow()) {
      auto &e = *graph->get_edge(v.get_edge_neighbors()[0]);

      auto &R = v.get_vessel_tree_data().resistances;

      if (!e.has_embedding_data())
        throw std::runtime_error("cannot determine coupling data for an unembedded graph");

      Point p = e.is_pointing_to(v_id) ? e.get_embedding_data().points.back() : e.get_embedding_data().points.front();

      // 1e3 since we have to convert kg -> g:
      auto p_out = pressure_values[v_id] * 1e3;

      // 1e3 since we have to convert kg -> g:
      auto R2 = R.back() * 1e3;

      auto radius_last = calculate_edge_tree_parameters(e).radii.back();
      auto radius_first = e.get_physical_data().radius;

      auto concentration = concentration_values[v_id];

      results.push_back({p, v.get_id(), p_out, concentration, R2, radius_last, radius_first, R.size()});
    }
  }

  return results;
}

void LinearizedHeartToBreast1DSolver::update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets) {
  for (auto v_id : graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &v = *graph->get_vertex(v_id);
    if (v.is_vessel_tree_outflow())
      // convert [Ba] to [kg / cm / s^2]
      v.update_vessel_tip_pressures(pressures_at_outlets.at(v_id) * 1e-3);
  }
}

void LinearizedHeartToBreast1DSolver::update_vessel_tip_concentrations(const std::map<size_t, double> &concentrations_at_outlets) {
  for (auto v_id : graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &v = *graph->get_vertex(v_id);
    if (v.is_vessel_tree_outflow())
      transport_solver->set_inflow_value(*graph, v, concentrations_at_outlets.at(v_id));
  }
}

void LinearizedHeartToBreast1DSolver::setup_graphs() {
  EmbeddedGraphReader graph_reader;

  graph_reader.append(path_linear_geometry, *graph);

  auto inflow_pressures = read_input_pressures(path_inflow_pressures);
  for (const auto &inflow_pressure : inflow_pressures) {
    auto &v_in = *graph->find_vertex_by_name(inflow_pressure.name);
    v_in.set_to_inflow_with_fixed_pressure(piecewise_linear_source_function(inflow_pressure.t, inflow_pressure.p, inflow_pressure.periodic));
  }

  set_0d_tree_boundary_conditions(graph);

  graph->finalize_bcs();

  flow_mesh_partitioner(d_comm, *graph, d_degree);
}

void LinearizedHeartToBreast1DSolver::setup_solver(size_t degree, double tau) {
  setup_solver_flow(degree, tau);
  setup_solver_transport(degree);
}

void LinearizedHeartToBreast1DSolver::setup_solver_flow(size_t degree, double tau) {
  d_degree = degree;
  d_tau_flow = tau;
  d_dof_map = std::make_shared<DofMap>(*graph);
  d_dof_map->create(d_comm, *graph, 2, d_degree, true);
  solver = std::make_shared<ImplicitLinearFlowSolver>(d_comm, graph, d_dof_map, d_degree);
  solver->setup(tau);
}

void LinearizedHeartToBreast1DSolver::setup_solver_transport(size_t degree) {
  // upwind evaluators:
  upwind_evaluator = std::make_shared<LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, d_dof_map);
  upwind_provider = std::make_shared<UpwindProviderLinearizedFlow>(graph, upwind_evaluator, solver);

  auto dof_map_transport = std::make_shared<DofMap>(*graph);
  DofMap::create_for_transport(MPI_COMM_WORLD, {graph}, {dof_map_transport}, degree);

  transport_solver = std::make_shared<ImplicitTransportSolver>(MPI_COMM_WORLD,
                                                               std::vector<std::shared_ptr<GraphStorage>>({graph}),
                                                               std::vector<std::shared_ptr<DofMap>>({dof_map_transport}),
                                                               std::vector<std::shared_ptr<UpwindProvider>>({upwind_provider}),
                                                               degree);
}


void LinearizedHeartToBreast1DSolver::setup_output() {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto dof_map_transport = transport_solver->get_dof_maps_transport().at(0);

  csv_writer = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv, graph);
  csv_writer->add_setup_data(d_dof_map, solver->p_component, "p");
  csv_writer->add_setup_data(d_dof_map, solver->q_component, "q");
  csv_writer->add_setup_data(dof_map_transport, 0, "c");
  csv_writer->setup();

  graph_pvd_writer = std::make_shared<GraphPVDWriter>(d_comm, output_folder_name, filename_pvd);

  vessel_tip_writer = std::make_shared<CSVVesselTipWriter>(MPI_COMM_WORLD,
                                                           output_folder_name, filename_csv_tips,
                                                           graph,
                                                           std::vector<std::shared_ptr<DofMap>>({d_dof_map, dof_map_transport, transport_solver->get_dof_maps_volume().back()}),
                                                           std::vector<std::string>({"p", "c", "V"}));

  // vessels ids and radii do not change, thus we can precalculate them
  fill_with_vessel_id(d_comm, *graph, points, vessel_ids);
  fill_with_radius(d_comm, *graph, points, vessel_radii);
  fill_with_vessel_A0(d_comm, *graph, points, vessel_A);
}

} // namespace macrocirculation