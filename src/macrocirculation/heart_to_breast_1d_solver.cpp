////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_1d_solver.hpp"

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
#include "interpolate_to_vertices.hpp"
#include "nonlinear_linear_coupling.hpp"
#include "tip_vertex_dof_integrator.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {


HeartToBreast1DSolver::HeartToBreast1DSolver(MPI_Comm comm)
    : d_comm(comm),
      d_degree(2),
      graph_nl{std::make_shared<GraphStorage>()},
      graph_li{std::make_shared<GraphStorage>()},
      coupling{std::make_shared<NonlinearLinearCoupling>(d_comm, graph_nl, graph_li)},
      d_tau(NAN),
      t(0),
      solver{nullptr},
      d_integrator_running(false) {}

void HeartToBreast1DSolver::start_0d_pressure_integrator() {
  integrator->reset();
  d_integrator_running = true;
}

std::vector<VesselTipCouplingData> HeartToBreast1DSolver::stop_0d_pressure_integrator() {
  d_integrator_running = false;

  auto values = integrator->get_integral_value({last_arterial_tip_index,
                                                first_vene_tip_index});

  std::vector<VesselTipCouplingData> results;

  for (auto v_id : graph_li->get_vertex_ids()) {
    auto &v = *graph_li->get_vertex(v_id);

    if (v.is_vessel_tree_outflow()) {
      auto &e = *graph_li->get_edge(v.get_edge_neighbors()[0]);

      if (!e.has_embedding_data())
        throw std::runtime_error("cannot determine coupling data for an unembedded graph");

      if (v.get_vessel_tree_data().resistances.size() != 7)
        throw std::runtime_error("unexpected number of resistors");

      Point p = e.is_pointing_to(v_id) ? e.get_embedding_data().points.back() : e.get_embedding_data().points.front();

      auto p_art = values[v_id][0];
      auto p_ven = values[v_id][1];

      results.push_back({p, p_art, p_ven});
    }
  }

  return results;
}

void HeartToBreast1DSolver::setup_graphs() {
  EmbeddedGraphReader graph_reader;

  graph_reader.append(path_nonlinear_geometry, *graph_nl);

  auto &v_in = *graph_nl->find_vertex_by_name("cw_in");
  v_in.set_to_inflow(heart_beat_inflow(485.0));

  graph_reader.append(path_linear_geometry, *graph_li);
  graph_reader.set_boundary_data(path_boundary_linear, *graph_li);
  graph_reader.set_boundary_data(path_boundary_nonlinear, *graph_nl);

  convert_rcr_to_partitioned_tree_bcs(graph_li);

  coupling->add_coupled_vertices("cw_out_1_1", "bg_132");
  coupling->add_coupled_vertices("cw_out_1_2", "bg_141");
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

void HeartToBreast1DSolver::setup_solver(size_t degree, double tau) {
  d_degree = degree;
  d_tau = tau;
  solver = std::make_shared<CoupledExplicitImplicit1DSolver>(d_comm, coupling, graph_nl, graph_li, d_degree, d_degree);
  solver->setup(tau);
  solver->get_explicit_solver()->use_ssp_method();

  auto dof_map_li = solver->get_implicit_dof_map();
  integrator = std::make_shared<TipVertexDofIntegrator>(d_comm, graph_li, dof_map_li);
}

void HeartToBreast1DSolver::setup_output() {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto dof_map_li = solver->get_implicit_dof_map();
  auto dof_map_nl = solver->get_explicit_dof_map();

  csv_writer_nl = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv_nl, graph_nl);
  csv_writer_nl->add_setup_data(dof_map_nl, get_solver_nl().A_component, "a");
  csv_writer_nl->add_setup_data(dof_map_nl, get_solver_nl().Q_component, "q");
  csv_writer_nl->setup();

  csv_writer_li = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv_li, graph_li);
  csv_writer_li->add_setup_data(dof_map_li, get_solver_li().p_component, "p");
  csv_writer_li->add_setup_data(dof_map_li, get_solver_li().q_component, "q");
  csv_writer_li->setup();

  graph_pvd_writer = std::make_shared<GraphPVDWriter>(d_comm, output_folder_name, filename_pvd);

  vessel_tip_writer = std::make_shared<CSVVesselTipWriter>(d_comm, output_folder_name, filename_csv_tips, graph_li, dof_map_li);

  // vessels ids do not change, thus we can precalculate them
  fill_with_vessel_id(d_comm, *graph_li, points, vessel_ids_li);
}

void HeartToBreast1DSolver::write_output() {
  if (solver == nullptr)
    throw std::runtime_error("solver must be initialized before output");

  auto dof_map_li = solver->get_implicit_dof_map();
  auto dof_map_nl = solver->get_explicit_dof_map();

  csv_writer_nl->add_data("a", get_solver_nl().get_solution());
  csv_writer_nl->add_data("q", get_solver_nl().get_solution());
  csv_writer_nl->write(t);

  csv_writer_li->add_data("p", get_solver_li().get_solution());
  csv_writer_li->add_data("q", get_solver_li().get_solution());
  csv_writer_li->write(t);

  interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, get_solver_li().p_component, get_solver_li().get_solution(), points, p_vertex_values);
  interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, get_solver_li().q_component, get_solver_li().get_solution(), points, q_vertex_values);

  graph_pvd_writer->set_points(points);
  graph_pvd_writer->add_vertex_data("p", p_vertex_values);
  graph_pvd_writer->add_vertex_data("q", q_vertex_values);
  graph_pvd_writer->add_vertex_data("vessel_id", vessel_ids_li);
  graph_pvd_writer->write(t);

  vessel_tip_writer->write(t, get_solver_li().get_solution());
}

void HeartToBreast1DSolver::setup(size_t degree, double tau) {
  setup_graphs();
  setup_solver(degree, tau);
  setup_output();
}

void HeartToBreast1DSolver::solve() {
  if (d_integrator_running)
    integrator->update_vertex_dof(get_solver_li().get_solution(), d_tau);
  solver->solve(d_tau, t);
  t += d_tau;
}

double HeartToBreast1DSolver::get_time() const { return t; }

} // namespace macrocirculation