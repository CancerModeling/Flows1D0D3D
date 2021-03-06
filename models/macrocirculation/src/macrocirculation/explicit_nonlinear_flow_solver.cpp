////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.hpp"
#include <cassert>

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"
#include "interpolate_to_vertices.hpp"
#include "quantities_of_interest.hpp"
#include "right_hand_side_evaluator.hpp"
#include "time_integrators.hpp"
#include "vessel_data_storage.hpp"
#include "communicator.hpp"
#include "communication/mpi.hpp"

namespace macrocirculation {

namespace lm = libMesh;

void interpolate_constant(const GraphStorage &graph, const DofMapNetwork &dof_map, double value, std::size_t component, std::vector<double> &result) {
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    dof_map.dof_indices(*edge, dof_indices, component);
    result[dof_indices[0]] = value;
  }
}

void set_to_A0(const GraphStorage &graph, const DofMapNetwork &dof_map, const VesselDataStorage &vessel_data, std::vector<double> &result) {
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    const auto &data = vessel_data.get_parameters(*edge);
    // set Q
    dof_map.dof_indices(*edge, dof_indices, 0);
    result[dof_indices[0]] = 0;
    // set A
    dof_map.dof_indices(*edge, dof_indices, 1);
    result[dof_indices[0]] = data.A0;
  }
}

template<std::size_t degree>
ExplicitNonlinearFlowSolver<degree>::ExplicitNonlinearFlowSolver(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_vessel_data(std::move(vessel_data)),
      d_dof_map(std::make_shared<SimpleDofMapNetwork>(2, degree + 1, d_graph->num_edges())),
      d_communicator(std::make_shared<Communicator>(comm, d_graph, d_dof_map)),
      d_right_hand_side_evaluator(std::make_shared<RightHandSideEvaluator<degree>>(d_comm, d_graph, d_vessel_data, d_dof_map)),
      d_time_integrator(std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof())),
      d_tau(2.5e-4 / 4),
      d_t_now(0),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof()) {
  // set A constant to A0
  set_to_A0(*d_graph, *d_dof_map, *d_vessel_data, d_u_prev);
  set_to_A0(*d_graph, *d_dof_map, *d_vessel_data, d_u_now);
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
double ExplicitNonlinearFlowSolver<degree>::get_time() const { return d_t_now; }

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::set_tau(double tau) { d_tau = tau; }

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::use_explicit_euler_method() {
  d_time_integrator = std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof());
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::use_ssp_method() {
  d_time_integrator = std::make_unique<TimeIntegrator>(create_ssp_method(), d_dof_map->num_dof());
}

template<std::size_t degree>
RightHandSideEvaluator<degree> &ExplicitNonlinearFlowSolver<degree>::get_rhs_evaluator() { return *d_right_hand_side_evaluator; }

template<std::size_t degree>
DofMapNetwork &ExplicitNonlinearFlowSolver<degree>::get_dof_map() { return *d_dof_map; }

template<std::size_t degree>
Communicator &ExplicitNonlinearFlowSolver<degree>::get_communicator() { return *d_communicator; }

template<std::size_t degree>
std::vector<double> &ExplicitNonlinearFlowSolver<degree>::get_solution() { return d_u_now; }

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_solution_on_vertices(std::vector< Point >& points, std::vector<double> &Q_values, std::vector<double> &A_values) const {
  assert(points.size() == d_graph->num_edges() * 2);
  assert(Q_values.size() == d_graph->num_edges() * 2);
  assert(A_values.size() == d_graph->num_edges() * 2);

  interpolate_to_vertices<degree>(d_comm, *d_graph, *d_dof_map, 0, d_u_now, points, Q_values);
  interpolate_to_vertices<degree>(d_comm, *d_graph, *d_dof_map, 1, d_u_now, points, A_values);
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_total_pressure_on_vertices(std::vector< Point >& points, std::vector<double> &p_values) const {
  assert(p_values.size() == d_graph->num_edges() * 2);

  calculate_total_pressure<degree>(d_comm, *d_graph, *d_vessel_data, *d_dof_map, d_u_now, points, p_values);
}

template<std::size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_static_pressure_on_vertices(std::vector< Point >& points, std::vector<double> &p_values) const {
  assert(p_values.size() == d_graph->num_edges() * 2);

  calculate_static_pressure<degree>(d_comm, *d_graph, *d_vessel_data, *d_dof_map, d_u_now, points, p_values);
}

// template instantiations:
template class ExplicitNonlinearFlowSolver<0>;
template class ExplicitNonlinearFlowSolver<1>;
template class ExplicitNonlinearFlowSolver<2>;
template class ExplicitNonlinearFlowSolver<3>;

} // namespace macrocirculation