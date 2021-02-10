////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.h"
#include <cassert>

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "gmm.h"
#include "graph_storage.hpp"
#include "interpolate_to_vertices.hpp"
#include "right_hand_side_evaluator.hpp"

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

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph)
    : d_graph(std::move(graph)),
      d_dof_map(std::make_shared<SimpleDofMapNetwork>(2, degree + 1, d_graph->num_edges())),
      d_right_hand_side_evaluator(std::make_shared<RightHandSideEvaluator<degree>>(d_graph, d_dof_map)),
      d_tau(2.5e-4 / 4),
      d_t_now(0),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof()),
      d_k0(d_dof_map->num_dof()),
      d_A0(6.97) // 6.97 cm^2,      TODO: Check if units are consistent!
{
  // set A constant to A0
  interpolate_constant(*d_graph, *d_dof_map, d_A0, 1, d_u_prev);
  interpolate_constant(*d_graph, *d_dof_map, d_A0, 1, d_u_now);
}

void ExplicitNonlinearFlowSolver::solve() {
  d_u_prev = d_u_now;
  d_t_now += d_tau;
  const double t_prev = d_t_now - d_tau;

  d_right_hand_side_evaluator->evaluate(t_prev, d_u_prev, d_k0);

  // u_now = u_prev + tau * k0
  gmm::add(d_u_prev, gmm::scaled(d_k0, d_tau), d_u_now);
}

double ExplicitNonlinearFlowSolver::get_time() const { return d_t_now; }

void ExplicitNonlinearFlowSolver::set_tau(double tau) { d_tau = tau; }

void ExplicitNonlinearFlowSolver::get_solution_on_vertices(std::vector<double> &Q_values, std::vector<double> &A_values) const {
  assert(Q_values.size() == d_graph->num_edges() * 2);
  assert(A_values.size() == d_graph->num_edges() * 2);

  FETypeNetwork<degree> fe(create_trapezoidal_rule());

  interpolate_to_vertices(*d_graph, *d_dof_map, fe, 0, d_u_now, Q_values);
  interpolate_to_vertices(*d_graph, *d_dof_map, fe, 1, d_u_now, A_values);
}

} // namespace macrocirculation