////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_1d_solver.hpp"

#include <utility>

#include "0d_boundary_conditions.hpp"
#include "communication/mpi.hpp"
#include "coupled_explicit_implicit_1d_solver.hpp"
#include "csv_vessel_tip_writer.hpp"
#include "dof_map.hpp"
#include "embedded_graph_reader.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "graph_csv_writer.hpp"
#include "graph_partitioner.hpp"
#include "graph_pvd_writer.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "implicit_transport_solver.hpp"
#include "interpolate_to_vertices.hpp"
#include "linearized_flow_upwind_evaluator.hpp"
#include "nonlinear_flow_upwind_evaluator.hpp"
#include "nonlinear_linear_coupling.hpp"
#include "tip_vertex_dof_integrator.hpp"
#include "vessel_formulas.hpp"
#include "vessel_tree_flow_integrator.hpp"

namespace macrocirculation {


HeartToBreast1DSolver::HeartToBreast1DSolver(MPI_Comm comm)
    : d_comm(comm),
      d_degree(2),
      d_is_setup(false),
      graph_nl{std::make_shared<GraphStorage>()},
      graph_li{std::make_shared<GraphStorage>()},
      coupling{std::make_shared<NonlinearLinearCoupling>(d_comm, graph_nl, graph_li)},
      d_tau_flow{NAN},
      solver{nullptr},
      d_integrator_running(false),
      d_flow_integrator_running(false) {}

void HeartToBreast1DSolver::start_0d_pressure_integrator() {
  integrator->reset();
  d_integrator_running = true;
}

void HeartToBreast1DSolver::start_0d_flow_integrator() {
  flow_integrator->reset();
  d_flow_integrator_running = true;
}

std::vector<VesselTreeFlowIntegratorResult> HeartToBreast1DSolver::stop_0d_flow_integrator() {
  d_flow_integrator_running = false;
  return flow_integrator->calculate();
}

void HeartToBreast1DSolver::stop_0d_flow_integrator_and_write() {
  d_flow_integrator_running = false;
  flow_integrator->write(output_folder_name, "0d_averaged_flows");
}

std::vector<VesselTipAverageCouplingData> HeartToBreast1DSolver::stop_0d_pressure_integrator() {
  d_integrator_running = false;

  auto values = integrator->get_integral_value({last_arterial_tip_index, first_vene_tip_index});

  double time_interval = integrator->get_integration_time();

  std::vector<VesselTipAverageCouplingData> results;

  for (auto v_id : graph_li->get_vertex_ids()) {
    auto &v = *graph_li->get_vertex(v_id);

    if (v.is_vessel_tree_outflow()) {
      auto &e = *graph_li->get_edge(v.get_edge_neighbors()[0]);

      auto &R = v.get_vessel_tree_data().resistances;

      if (!e.has_embedding_data())
        throw std::runtime_error("cannot determine coupling data for an unembedded graph");

      Point p = e.is_pointing_to(v_id) ? e.get_embedding_data().points.back() : e.get_embedding_data().points.front();

      // 1e3 since we have to convert kg -> g:
      auto p_art = values[v_id][0] * 1e3 / time_interval;
      auto p_ven = values[v_id][1] * 1e3 / time_interval;

      // 1e3 since we have to convert kg -> g:
      // auto R1 = calculate_R1(e.get_physical_data());
      // auto R2_art = (R[last_arterial_tip_index] - R1) * 1e3;
      // auto R2_cap = (R[capillary_tip_index] - R1) * 1e3;

      // 1e3 since we have to convert kg -> g:
      auto R2_art = (R[last_arterial_tip_index]) * 1e3;
      auto R2_cap = (R[capillary_tip_index]) * 1e3;

      results.push_back({p, v.get_id(), p_art, p_ven, R2_art, R2_cap});
    }
  }

  return results;
}

std::map<size_t, double> get_vessel_tip_dof_values(MPI_Comm comm,
                                                   const GraphStorage &graph,
                                                   const DofMap &dof_map,
                                                   const PetscVec &u) {
  std::map<size_t, double> data;
  std::vector<double> vertex_dof_values;
  for (auto v_id : graph.get_vertex_ids()) {
    auto v = graph.get_vertex(v_id);
    if (v->is_vessel_tree_outflow()) {
      data[v_id] = 0;
      auto &e = *graph.get_edge(v->get_edge_neighbors()[0]);
      if (e.rank() == mpi::rank(comm)) {
        auto local_dof_map = dof_map.get_local_dof_map(*v);
        const auto &dof_indices = local_dof_map.dof_indices();
        vertex_dof_values.resize(dof_indices.size());
        extract_dof(dof_indices, u, vertex_dof_values);
        data.at(v_id) = vertex_dof_values.back();
      }
      MPI_Bcast(&data[v_id], 1, MPI_DOUBLE, e.rank(), comm);
    }
  }
  return data;
}

void HeartToBreast1DSolver::set_path_inflow_pressures(const std::string &path){
  if (d_is_setup)
    throw std::runtime_error("HeartToBreast1DSolver::set_path_inflow_pressures: cannot be called after setup.");

  path_inflow_pressures = path;
}

void HeartToBreast1DSolver::set_path_nonlinear_geometry(const std::string& path) {
  if (d_is_setup)
    throw std::runtime_error("HeartToBreast1DSolver::set_path_nonlinear_geometry: cannot be called after setup.");

  path_nonlinear_geometry = path;
}

std::vector<VesselTipCurrentCouplingData> HeartToBreast1DSolver::get_vessel_tip_pressures() {
  auto &dof_map_flow = *solver->get_implicit_dof_map();
  auto &u_flow = solver->get_implicit_solver()->get_solution();
  auto pressure_values = get_vessel_tip_dof_values(d_comm, *graph_li, dof_map_flow, u_flow);

  auto &dof_map_transport = *transport_solver->get_dof_maps_transport().back();
  auto &u_transport = transport_solver->get_solution();
  auto concentration_values = get_vessel_tip_dof_values(d_comm, *graph_li, dof_map_transport, u_transport);

  std::vector<VesselTipCurrentCouplingData> results;

  for (auto v_id : graph_li->get_vertex_ids()) {
    auto &v = *graph_li->get_vertex(v_id);

    if (v.is_vessel_tree_outflow()) {
      auto &e = *graph_li->get_edge(v.get_edge_neighbors()[0]);

      auto &R = v.get_vessel_tree_data().resistances;

      if (!e.has_embedding_data())
        throw std::runtime_error("cannot determine coupling data for an unembedded graph");

      Point p = e.is_pointing_to(v_id) ? e.get_embedding_data().points.back() : e.get_embedding_data().points.front();

      // 1e3 since we have to convert kg -> g:
      auto p_out = pressure_values[v_id] * 1e3;

      // 1e3 since we have to convert kg -> g:
      auto R2 = R.back() * 1e3;

      auto radius = calculate_edge_tree_parameters(e).radii.back();

      auto concentration = concentration_values[v_id];

      results.push_back({p, v.get_id(), p_out, concentration, R2, radius, R.size()});
    }
  }

  return results;
}

void HeartToBreast1DSolver::update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets) {
  for (auto v_id : graph_li->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &v = *graph_li->get_vertex(v_id);
    if (v.is_vessel_tree_outflow())
      // convert [Ba] to [kg / cm / s^2]
      v.update_vessel_tip_pressures(pressures_at_outlets.at(v_id) * 1e-3);
  }
}

void HeartToBreast1DSolver::update_vessel_tip_concentrations(const std::map<size_t, double> &concentrations_at_outlets) {
  for (auto v_id : graph_li->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &v = *graph_li->get_vertex(v_id);
    if (v.is_vessel_tree_outflow())
      transport_solver->set_inflow_value(*graph_li, v, concentrations_at_outlets.at(v_id));
  }
}

void HeartToBreast1DSolver::setup_graphs(BoundaryModel bmodel) {
  EmbeddedGraphReader graph_reader;

  graph_reader.append(path_nonlinear_geometry, *graph_nl);

  graph_reader.append(path_nonlinear_geometry, *graph_nl);
  if (!path_inflow_pressures.empty()) {
    auto inflow_pressures = read_input_pressures(path_inflow_pressures);
    for (const auto& inflow_pressure : inflow_pressures)
    {
      auto &v_in = *graph_nl->find_vertex_by_name(inflow_pressure.name);
      v_in.set_to_inflow_with_fixed_pressure(piecewise_linear_source_function(inflow_pressure.t, inflow_pressure.p));
    }
  } else {
    auto &v_in = *graph_nl->find_vertex_by_name("cw_in");
    v_in.set_to_inflow_with_fixed_flow(heart_beat_inflow(485.0));
    graph_reader.set_boundary_data(path_boundary_nonlinear, *graph_nl);
  };

  graph_reader.append(path_linear_geometry, *graph_li);
  graph_reader.set_boundary_data(path_boundary_linear, *graph_li);

  if (bmodel == BoundaryModel::DiscreteRCRChain)
    convert_rcr_to_partitioned_tree_bcs(graph_li);
  else if (bmodel == BoundaryModel::DiscreteRCRTree)
    set_0d_tree_boundary_conditions(graph_li, "bg");
  else
    throw std::runtime_error("Boundary model type was not implemented yet.");

  coupling->add_coupled_vertices("cw_out_1_1", "bg_141");
  coupling->add_coupled_vertices("cw_out_1_2", "bg_139");
  coupling->add_coupled_vertices("cw_out_1_3", "bg_132");
  coupling->add_coupled_vertices("cw_out_2_1", "bg_135");
  coupling->add_coupled_vertices("cw_out_2_2", "bg_119");

  graph_nl->finalize_bcs();
  graph_li->finalize_bcs();

  // we assume a degree of 2 for partitioning
  flow_mesh_partitioner(d_comm, *graph_nl, d_degree);
  flow_mesh_partitioner(d_comm, *graph_li, d_degree);
}

ImplicitLinearFlowSolver &HeartToBreast1DSolver::get_solver_li() {
  if (solver == nullptr)
    throw std::runtime_error("solver not initialized");
  return *solver->get_implicit_solver();
}

ExplicitNonlinearFlowSolver &HeartToBreast1DSolver::get_solver_nl() {
  if (solver == nullptr)
    throw std::runtime_error("solver not initialized");
  return *solver->get_explicit_solver();
}

void HeartToBreast1DSolver::setup_solver_flow(size_t degree, double tau_flow) {
  d_degree = degree;
  d_tau_flow = tau_flow;
  solver = std::make_shared<CoupledExplicitImplicit1DSolver>(d_comm, coupling, graph_nl, graph_li, d_degree, d_degree);
  solver->setup(tau_flow);
  solver->get_explicit_solver()->use_ssp_method();

  auto dof_map_li = solver->get_implicit_dof_map();
  integrator = std::make_shared<TipVertexDofIntegrator>(d_comm, graph_li, dof_map_li);
  flow_integrator = std::make_shared<VesselTreeFlowIntegrator>(d_comm, graph_li, dof_map_li);
}

void HeartToBreast1DSolver::setup_solver(size_t degree, double tau) {
  setup_solver_flow(degree, tau);
  setup_solver_transport(degree);
}

void HeartToBreast1DSolver::setup_solver_transport(size_t degree) {
  auto dof_map_li = solver->get_implicit_dof_map();
  auto dof_map_nl = solver->get_explicit_dof_map();

  // upwind evaluators:
  upwind_evaluator_nl = std::make_shared<NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_nl, dof_map_nl);
  upwind_provider_nl = std::make_shared<UpwindProviderNonlinearFlow>(upwind_evaluator_nl, solver->get_explicit_solver());

  // upwind evaluators:
  upwind_evaluator_li = std::make_shared<LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_li, dof_map_li);
  upwind_provider_li = std::make_shared<UpwindProviderLinearizedFlow>(graph_li, upwind_evaluator_li, solver->get_implicit_solver());

  auto dof_map_transport_nl = std::make_shared<DofMap>(*graph_nl);
  auto dof_map_transport_li = std::make_shared<DofMap>(*graph_li);
  DofMap::create_for_transport(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, degree);

  transport_solver = std::make_shared<ImplicitTransportSolver>(MPI_COMM_WORLD,
                                                               std::vector<std::shared_ptr<GraphStorage>>({graph_nl, graph_li}),
                                                               std::vector<std::shared_ptr<DofMap>>({dof_map_transport_nl, dof_map_transport_li}),
                                                               std::vector<std::shared_ptr<UpwindProvider>>({upwind_provider_nl, upwind_provider_li}),
                                                               degree);
}

void HeartToBreast1DSolver::setup_output() {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto dof_map_li = solver->get_implicit_dof_map();
  auto dof_map_nl = solver->get_explicit_dof_map();

  auto dof_map_transport_nl = transport_solver->get_dof_maps_transport().at(0);
  auto dof_map_transport_li = transport_solver->get_dof_maps_transport().at(1);

  csv_writer_nl = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv_nl, graph_nl);
  csv_writer_nl->add_setup_data(dof_map_nl, get_solver_nl().A_component, "a");
  csv_writer_nl->add_setup_data(dof_map_nl, get_solver_nl().Q_component, "q");
  csv_writer_nl->add_setup_data(dof_map_transport_nl, 0, "c");
  csv_writer_nl->setup();

  csv_writer_li = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv_li, graph_li);
  csv_writer_li->add_setup_data(dof_map_li, get_solver_li().p_component, "p");
  csv_writer_li->add_setup_data(dof_map_li, get_solver_li().q_component, "q");
  csv_writer_li->add_setup_data(dof_map_transport_li, 0, "c");
  csv_writer_li->setup();

  graph_pvd_writer = std::make_shared<GraphPVDWriter>(d_comm, output_folder_name, filename_pvd);

  vessel_tip_writer_nl = std::make_shared<CSVVesselTipWriter>(MPI_COMM_WORLD,
                                                              output_folder_name, filename_csv_tips_nl,
                                                              graph_nl,
                                                              std::vector<std::shared_ptr<DofMap>>({dof_map_nl, dof_map_transport_nl, transport_solver->get_dof_maps_volume().front()}),
                                                              std::vector<std::string>({"p", "c", "V"}));

  vessel_tip_writer_li = std::make_shared<CSVVesselTipWriter>(MPI_COMM_WORLD,
                                                              output_folder_name, filename_csv_tips_li,
                                                              graph_li,
                                                              std::vector<std::shared_ptr<DofMap>>({dof_map_li, dof_map_transport_li, transport_solver->get_dof_maps_volume().back()}),
                                                              std::vector<std::string>({"p", "c", "V"}));

  // vessels ids and radii do not change, thus we can precalculate them
  fill_with_vessel_id(d_comm, *graph_li, points, vessel_ids_li);
  fill_with_radius(d_comm, *graph_li, points, vessel_radii_li);
  fill_with_vessel_A0(d_comm, *graph_li, points, vessel_A_li);
}

void HeartToBreast1DSolver::write_output(double t) {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto dof_map_li = solver->get_implicit_dof_map();
  auto dof_map_nl = solver->get_explicit_dof_map();

  auto dof_map_transport_li = transport_solver->get_dof_maps_transport().at(1);

  csv_writer_nl->add_data("a", get_solver_nl().get_solution());
  csv_writer_nl->add_data("q", get_solver_nl().get_solution());
  csv_writer_nl->add_data("c", transport_solver->get_solution());
  csv_writer_nl->write(t);

  csv_writer_li->add_data("p", get_solver_li().get_solution());
  csv_writer_li->add_data("q", get_solver_li().get_solution());
  csv_writer_li->add_data("c", transport_solver->get_solution());
  csv_writer_li->write(t);

  interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, get_solver_li().p_component, get_solver_li().get_solution(), points, p_vertex_values);
  interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, get_solver_li().q_component, get_solver_li().get_solution(), points, q_vertex_values);
  interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_transport_li, 0, transport_solver->get_solution(), points, c_vertex_values);

  graph_pvd_writer->set_points(points);
  graph_pvd_writer->add_vertex_data("p", p_vertex_values);
  graph_pvd_writer->add_vertex_data("q", q_vertex_values);
  graph_pvd_writer->add_vertex_data("c", c_vertex_values);
  graph_pvd_writer->add_vertex_data("vessel_id", vessel_ids_li);
  graph_pvd_writer->add_vertex_data("r", vessel_radii_li);
  graph_pvd_writer->add_vertex_data("A", vessel_A_li);
  graph_pvd_writer->write(t);

  vessel_tip_writer_nl->write(t, {get_solver_nl().get_solution()}, {transport_solver->get_solution(), transport_solver->get_volumes()});
  vessel_tip_writer_li->write(t, {get_solver_li().get_solution(), transport_solver->get_solution(), transport_solver->get_volumes()});
}

void HeartToBreast1DSolver::setup(size_t degree, double tau, BoundaryModel boundary_model) {
  d_is_setup = true;
  setup_graphs(boundary_model);
  setup_solver(degree, tau);
  setup_output();
}

void HeartToBreast1DSolver::solve_flow(double tau, double t) {
  if (std::abs(tau - d_tau_flow) > 1e-12)
    throw std::runtime_error("changing time step width not supported yet.");

  if (d_integrator_running)
    integrator->update_vertex_dof(get_solver_li().get_solution(), tau);
  if (d_flow_integrator_running)
    flow_integrator->add(get_solver_li().get_solution(), tau);
  solver->solve(tau, t);
}

void HeartToBreast1DSolver::solve_transport(double tau, double t) {
  upwind_evaluator_nl->init(t, solver->get_explicit_solver()->get_solution());
  upwind_evaluator_li->init(t, solver->get_implicit_solver()->get_solution());
  transport_solver->solve(tau, t);
}

void HeartToBreast1DSolver::apply_slope_limiter_transport(double t) {
  transport_solver->apply_slope_limiter(t);
}

void HeartToBreast1DSolver::set_output_folder(std::string output_dir) { output_folder_name = std::move(output_dir); }

} // namespace macrocirculation