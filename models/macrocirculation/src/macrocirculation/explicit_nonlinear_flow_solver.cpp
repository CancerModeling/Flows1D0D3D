////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.hpp"

#include <cassert>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "graph_storage.hpp"
#include "right_hand_side_evaluator.hpp"
#include "time_integrators.hpp"

namespace macrocirculation {

void interpolate_constant(MPI_Comm comm,
                          const GraphStorage &graph,
                          const DofMap &dof_map,
                          double value,
                          std::size_t component,
                          std::vector<double> &result) {
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);
    const auto local_dof_map = dof_map.get_local_dof_map(*edge);
    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
    }
    result[dof_indices[0]] = value;
  }
}

void set_to_A0(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, std::vector<double> &result) {
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);
    assert(edge->has_physical_data());
    const auto &data = edge->get_physical_data();
    const auto &local_dof_map = dof_map.get_local_dof_map(*edge);
    dof_indices.resize(local_dof_map.num_basis_functions());
    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      // set Q
      local_dof_map.dof_indices(micro_edge_id, 0, dof_indices);
      result[dof_indices[0]] = 0;
      // set A
      local_dof_map.dof_indices(micro_edge_id, 1, dof_indices);
      result[dof_indices[0]] = data.A0;
    }
  }
}

template<std::size_t degree>
ExplicitNonlinearFlowSolver<degree>::ExplicitNonlinearFlowSolver(MPI_Comm comm,
                                                                 std::shared_ptr<GraphStorage> graph,
                                                                 std::shared_ptr<DofMap> dof_map)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_right_hand_side_evaluator(std::make_shared<RightHandSideEvaluator>(d_comm, d_graph, d_dof_map, degree)),
      d_time_integrator(std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof())),
      d_tau(2.5e-4 / 4),
      d_t_now(0),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof()) {
  // set A constant to A0
  set_to_A0(d_comm, *d_graph, *d_dof_map, d_u_prev);
  set_to_A0(d_comm, *d_graph, *d_dof_map, d_u_now);
}

// we need the destructor here, to use unique_ptrs with forward declared classes.
template<std::size_t degree>
ExplicitNonlinearFlowSolver<degree>::~ExplicitNonlinearFlowSolver() = default;

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::solve() {
  d_u_prev = d_u_now;
  d_t_now += d_tau;
  const double t_prev = d_t_now - d_tau;
  d_time_integrator->apply<degree>(d_u_prev, t_prev, d_tau, *d_right_hand_side_evaluator, d_u_now);
}

template<std::size_t degree>
double ExplicitNonlinearFlowSolver<degree>::get_time() const {
  return d_t_now;
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::set_tau(double tau) {
  d_tau = tau;
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::use_explicit_euler_method() {
  d_time_integrator = std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof());
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::use_ssp_method() {
  d_time_integrator = std::make_unique<TimeIntegrator>(create_ssp_method(), d_dof_map->num_dof());
}

template<std::size_t degree>
RightHandSideEvaluator &ExplicitNonlinearFlowSolver<degree>::get_rhs_evaluator() {
  return *d_right_hand_side_evaluator;
}

template<std::size_t degree>
DofMap &ExplicitNonlinearFlowSolver<degree>::get_dof_map() {
  return *d_dof_map;
}

template<std::size_t degree>
std::vector<double> &ExplicitNonlinearFlowSolver<degree>::get_solution() {
  return d_u_now;
}

/*
template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_solution_on_vertices(std::vector<double> &Q_values, std::vector<double> &A_values) const {
  assert(Q_values.size() == d_graph->num_edges() * 2);
  assert(A_values.size() == d_graph->num_edges() * 2);

  FEType fe(create_trapezoidal_rule());

  // interpolate_to_vertices(*d_graph, *d_dof_map, fe, 0, d_u_now, Q_values);
  // interpolate_to_vertices(*d_graph, *d_dof_map, fe, 1, d_u_now, A_values);
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_total_pressure_on_vertices(std::vector<double> &p_values) const {
  assert(p_values.size() == d_graph->num_edges() * 2);

  FEType fe(create_trapezoidal_rule());

  calculate_total_pressure(*d_graph, *d_vessel_data, *d_dof_map, fe, d_u_now, p_values);
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_static_pressure_on_vertices(std::vector<double> &p_values) const {
  assert(p_values.size() == d_graph->num_edges() * 2);

  FETypeNetwork<degree> fe(create_trapezoidal_rule());

  calculate_total_pressure(*d_graph, *d_vessel_data, *d_dof_map, fe, d_u_now, p_values);
}
*/

// template instantiations:
template class ExplicitNonlinearFlowSolver<0>;
template class ExplicitNonlinearFlowSolver<1>;
template class ExplicitNonlinearFlowSolver<2>;
template class ExplicitNonlinearFlowSolver<3>;

} // namespace macrocirculation