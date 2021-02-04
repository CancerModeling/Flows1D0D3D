////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.h"

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

constexpr std::size_t degree = 2;

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph)
    : d_graph(std::move(graph)) {
  // TODO: Initialize the rest
}

void ExplicitNonlinearFlowSolver::solve() {
}

void ExplicitNonlinearFlowSolver::calculate_fluxes() {
  // dof indices for left and right edge
  std::vector< std::size_t > Q_dof_indices_l;
  std::vector< std::size_t > A_dof_indices_l;
  std::vector< std::size_t > Q_dof_indices_r;
  std::vector< std::size_t > A_dof_indices_r;

  // finite-element for left and right edge
  FETypeNetwork<degree> fe_l(create_trapezoidal_rule());
  FETypeNetwork<degree> fe_r(create_trapezoidal_rule());

  const std::size_t num_basis_functions = degree + 1;

  // local views of our previous solution
  std::vector< double > Q_prev_loc_r (num_basis_functions, 0);
  std::vector< double > A_prev_loc_r (num_basis_functions, 0);
  std::vector< double > Q_prev_loc_l (num_basis_functions, 0);
  std::vector< double > A_prev_loc_l (num_basis_functions, 0);

  // previous solution evaluated at quadrature points
  std::vector< double > Q_prev_qp_r (num_basis_functions, 0);
  std::vector< double > A_prev_qp_r (num_basis_functions, 0);
  std::vector< double > Q_prev_qp_l (num_basis_functions, 0);
  std::vector< double > A_prev_qp_l (num_basis_functions, 0);

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

      const double W2_l = get_W2_value(Q_l, A_l, d_G0, d_rho, d_A0);
      const double W1_r = get_W1_value(Q_r, A_r, d_G0, d_rho, d_A0);

    }
  }
}

} // namespace macrocirculation