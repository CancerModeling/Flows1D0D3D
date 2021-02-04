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
#include "graph_storage.hpp"

namespace macrocirculation {

constexpr std::size_t degree = 2;

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph)
    : d_graph(std::move(graph)),
      d_dof_map(std::make_shared<SimpleDofMapNetwork>(2, degree + 1, d_graph->num_edges())),
      d_tau(1e-2),
      d_t_now(0),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof()),
      d_Q_up(d_graph->num_vertices()),
      d_A_up(d_graph->num_vertices()),
      d_G0(592.4e2), // 592.4 10^2 Pa,  TODO: Check if units are consistent!
      d_rho(1.028),  // 1.028 g/cm^3,   TODO: Check if units are consistent!
      d_A0(6.97)     // 6.97 cm^2,      TODO: Check if units are consistent!
{
// TODO: Initialize the rest
}

void ExplicitNonlinearFlowSolver::solve() {
}

void ExplicitNonlinearFlowSolver::calculate_fluxes() {
  const std::size_t num_basis_functions = degree + 1;

  // dof indices for left and right edge
  std::vector< std::size_t > Q_dof_indices_l(num_basis_functions, 0);
  std::vector< std::size_t > A_dof_indices_l(num_basis_functions, 0);
  std::vector< std::size_t > Q_dof_indices_r(num_basis_functions, 0);
  std::vector< std::size_t > A_dof_indices_r(num_basis_functions, 0);

  // finite-element for left and right edge
  FETypeNetwork<degree> fe_l(create_trapezoidal_rule());
  FETypeNetwork<degree> fe_r(create_trapezoidal_rule());

  // local views of our previous solution
  std::vector< double > Q_prev_loc_r (num_basis_functions, 0);
  std::vector< double > A_prev_loc_r (num_basis_functions, 0);
  std::vector< double > Q_prev_loc_l (num_basis_functions, 0);
  std::vector< double > A_prev_loc_l (num_basis_functions, 0);

  // previous solution evaluated at quadrature points
  // we have 2 quadrature points for the trapezoidal rule
  std::vector< double > Q_prev_qp_r (2, 0);
  std::vector< double > A_prev_qp_r (2, 0);
  std::vector< double > Q_prev_qp_l (2, 0);
  std::vector< double > A_prev_qp_l (2, 0);

  for (const auto &v_id : d_graph->get_vertex_ids()) {
    const auto vertex = d_graph->get_vertex(v_id);

    // exterior boundary
    if (vertex->is_leaf()) {
      // TODO
    }
    // inner boundary
    else {
      const auto edge_l = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[1]);

      // we assert that the edge direction fits our assumptions
      assert(edge_l->get_vertex_neighbors()[1] == vertex->get_id());
      assert(edge_r->get_vertex_neighbors()[0] == vertex->get_id());

      fe_l.reinit(*edge_l);
      fe_r.reinit(*edge_r);

      d_dof_map->dof_indices(*edge_l, Q_dof_indices_l, 0);
      d_dof_map->dof_indices(*edge_l, A_dof_indices_l, 1);
      d_dof_map->dof_indices(*edge_r, Q_dof_indices_r, 0);
      d_dof_map->dof_indices(*edge_r, A_dof_indices_r, 1);

      extract_dof(Q_dof_indices_l, d_u_prev, Q_prev_loc_l);
      extract_dof(A_dof_indices_l, d_u_prev, A_prev_loc_l);
      extract_dof(Q_dof_indices_r, d_u_prev, Q_prev_loc_r);
      extract_dof(A_dof_indices_r, d_u_prev, A_prev_loc_r);

      fe_l.evaluate_dof_at_quadrature_points(Q_prev_loc_l,  Q_prev_qp_l);
      fe_l.evaluate_dof_at_quadrature_points(A_prev_loc_l,  A_prev_qp_l);
      fe_l.evaluate_dof_at_quadrature_points(Q_prev_loc_r,  Q_prev_qp_r);
      fe_l.evaluate_dof_at_quadrature_points(A_prev_loc_r,  A_prev_qp_r);

      const double Q_l = Q_prev_qp_l[1];
      const double A_l = A_prev_qp_l[1];
      const double Q_r = Q_prev_qp_r[0];
      const double A_r = A_prev_qp_r[0];

      const double W2_l = calculate_W2_value(Q_l, A_l, d_G0, d_rho, d_A0);
      const double W1_r = calculate_W1_value(Q_r, A_r, d_G0, d_rho, d_A0);

      double Q_up = 0, A_up = 0;
      solve_W12(Q_up, A_up, W1_r, W2_l, d_G0, d_rho, d_A0);

      d_Q_up[vertex->get_id()] = Q_up;
      d_A_up[vertex->get_id()] = A_up;
    }
  }
}

} // namespace macrocirculation