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
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "right_hand_side_evaluator.hpp"
#include "time_integrators.hpp"
#include "vessel_formulas.hpp"

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
    dof_indices.resize(local_dof_map.num_basis_functions());
    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
      result[dof_indices[0]] = value;
    }
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

  for (const auto &v_id : graph.get_active_vertex_ids(mpi::rank(comm))) {
    const auto vertex = graph.get_vertex(v_id);

    if (!vertex->is_windkessel_outflow())
      continue;

    assert(vertex->is_leaf());

    const auto edge = graph.get_edge(vertex->get_edge_neighbors()[0]);

    const auto &vertex_dof_map = dof_map.get_local_dof_map(*vertex);

    assert(vertex_dof_map.num_local_dof() == 1);

    const auto &data = edge->get_physical_data();
    const auto &vertex_dof_indices = vertex_dof_map.dof_indices();

    // set p
    result[vertex_dof_indices[0]] = calculate_p_from_QA(0, data.A0, data.G0, data.rho, data.A0);
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

template<size_t degree>
double ExplicitNonlinearFlowSolver<degree>::get_flow_at_vessel_tip(const Vertex &v) const {
  if (!v.is_leaf())
    throw std::runtime_error("flow can only be calculated at leafs");

  const auto e_id = v.get_edge_neighbors()[0];

  const auto &edge = *d_graph->get_edge(e_id);

  const auto ldofmap = d_dof_map->get_local_dof_map(edge);

  QuadratureFormula qf = create_gauss4();
  FETypeNetwork fe(qf, ldofmap.num_basis_functions() - 1);

  std::vector<size_t> dof_indices(ldofmap.num_basis_functions(), 0);

  if (edge.is_pointing_to(v.get_id())) {
    ldofmap.dof_indices(edge.num_micro_edges() - 1, 0, dof_indices);
  } else {
    ldofmap.dof_indices(0, 0, dof_indices);
  }

  std::vector<double> dof_values(ldofmap.num_basis_functions(), 0);
  extract_dof(dof_indices, d_u_now, dof_values);
  auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_values);

  if (edge.is_pointing_to(v.get_id())) {
    return boundary_values.right;
  } else {
    return -boundary_values.left;
  }
}

template<size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_1d_AQ_values_at_vertex(const Vertex &v, double &A, double &Q) const {
  if (!v.is_leaf())
    throw std::runtime_error("flow can only be calculated at leafs");

  const auto e_id = v.get_edge_neighbors()[0];

  const auto &edge = *d_graph->get_edge(e_id);

  const auto ldofmap = d_dof_map->get_local_dof_map(edge);

  QuadratureFormula qf = create_gauss4();
  FETypeNetwork fe(qf, ldofmap.num_basis_functions() - 1);

  std::vector<size_t> dof_indices_q(ldofmap.num_basis_functions(), 0);
  std::vector<size_t> dof_indices_a(ldofmap.num_basis_functions(), 0);

  if (edge.is_pointing_to(v.get_id())) {
    ldofmap.dof_indices(edge.num_micro_edges() - 1, 0, dof_indices_q);
    ldofmap.dof_indices(edge.num_micro_edges() - 1, 1, dof_indices_a);
  } else {
    ldofmap.dof_indices(0, 0, dof_indices_q);
    ldofmap.dof_indices(0, 1, dof_indices_a);
  }

  std::vector<double> dof_values_q(ldofmap.num_basis_functions(), 0);
  extract_dof(dof_indices_q, d_u_now, dof_values_q);
  auto boundary_values_q = fe.evaluate_dof_at_boundary_points(dof_values_q);

  std::vector<double> dof_values_a(ldofmap.num_basis_functions(), 0);
  extract_dof(dof_indices_a, d_u_now, dof_values_a);
  auto boundary_values_a = fe.evaluate_dof_at_boundary_points(dof_values_a);

  if (edge.is_pointing_to(v.get_id())) {
    Q = boundary_values_q.right;
    A = boundary_values_a.right;
  } else {
    Q = boundary_values_q.left;
    A = boundary_values_a.left;
  }
}

template<size_t degree>
void ExplicitNonlinearFlowSolver<degree>::get_1d_pq_values_at_vertex(const Vertex &v, double &p, double &q) const {
  auto &data = d_graph->get_edge(v.get_edge_neighbors()[0])->get_physical_data();

  double A, Q;
  get_1d_AQ_values_at_vertex(v, A, Q);

  q = Q;
  p = nonlinear::get_p_from_A(data, A);
}

template<size_t degree>
[[nodiscard]] Values0DModel ExplicitNonlinearFlowSolver<degree>::get_0D_values(const Vertex &v) const {
  Values0DModel result{0, 0};

  const auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[0]);

  if (edge.rank() == mpi::rank(MPI_COMM_WORLD)) {
    const auto &vertex_dof_map = d_dof_map->get_local_dof_map(v);
    const auto &vertex_dofs = vertex_dof_map.dof_indices();
    const auto p_c = d_u_now[vertex_dofs[0]];

    const auto &param = edge.get_physical_data();

    // TODO: Move this calculation to the vertex.
    const double R1 = param.rho * param.get_c0() / param.A0;
    const double R2 = v.get_peripheral_vessel_data().resistance - R1;

    const double q = (p_c - v.get_peripheral_vessel_data().p_out) / R2;

    result.p_c = p_c;
    result.q = q;
  }

  MPI_Bcast(&result.p_c, 1, MPI_DOUBLE, edge.rank(), MPI_COMM_WORLD);
  MPI_Bcast(&result.q, 1, MPI_DOUBLE, edge.rank(), MPI_COMM_WORLD);

  return result;
}

template<size_t degree>
void ExplicitNonlinearFlowSolver<degree>::evaluate_1d_AQ_values(const Edge &e, double s, double &A, double &Q) const {
  // on which micro edge is the given value
  auto micro_edge_id = static_cast<size_t>(std::ceil(e.num_micro_edges() * s));
  micro_edge_id = std::min(micro_edge_id, e.num_micro_edges() - 1);
  const double h = 1. / e.num_micro_edges();
  // get parametrization value on the micro edge
  const double s_tilde = 2 * (s - h * static_cast<double>(micro_edge_id)) / h - 1;
  // get dof information
  const auto &local_dof_map = d_dof_map->get_local_dof_map(e);
  std::vector<size_t> dof(local_dof_map.num_basis_functions());
  std::vector<double> dof_values(local_dof_map.num_basis_functions());
  local_dof_map.dof_indices(micro_edge_id, A_component, dof);
  extract_dof(dof, d_u_now, dof_values);
  A = FETypeNetwork::evaluate_dof(dof_values, s_tilde);
  local_dof_map.dof_indices(micro_edge_id, Q_component, dof);
  extract_dof(dof, d_u_now, dof_values);
  Q = FETypeNetwork::evaluate_dof(dof_values, s_tilde);
}

// template instantiations:
template class ExplicitNonlinearFlowSolver<0>;
template class ExplicitNonlinearFlowSolver<1>;
template class ExplicitNonlinearFlowSolver<2>;
template class ExplicitNonlinearFlowSolver<3>;

} // namespace macrocirculation