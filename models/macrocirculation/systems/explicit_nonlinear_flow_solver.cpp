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
#include "libmesh/libmesh.h"
#include "gmm.h"

namespace macrocirculation {

namespace lm = libMesh;

constexpr std::size_t degree = 2;

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph)
    : d_graph(std::move(graph)),
      d_dof_map(std::make_shared<SimpleDofMapNetwork>(2, degree + 1, d_graph->num_edges())),
      d_tau(1e-2),
      d_t_now(0),
      d_inflow_value_function(heart_beat_inflow),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof()),
      d_Q_up(d_graph->num_vertices()),
      d_A_up(d_graph->num_vertices()),
      d_rhs(d_dof_map->num_dof()),
      d_mass(d_dof_map->num_dof()),
      d_G0(592.4e2), // 592.4 10^2 Pa,  TODO: Check if units are consistent!
      d_rho(1.028),  // 1.028 kg/cm^3,  TODO: Check if units are consistent!
      d_A0(6.97),    // 6.97 cm^2,      TODO: Check if units are consistent!
      d_mu(4.5),     // 4.5 m Pa/s      TODO: Check if units are consistent!
      d_gamma(2)     // Poiseuille flow
{
  // set A constant to A0
  std::vector< std::size_t > dof_indices;
  for (const auto &e_id : d_graph->get_edge_ids()) {
    const auto edge = d_graph->get_edge(e_id);
    d_dof_map->dof_indices(*edge, dof_indices, 1);
    d_u_prev[dof_indices[0]] = d_A0;
    d_u_now[dof_indices[0]] = d_A0;
  }
}

void ExplicitNonlinearFlowSolver::solve() {
  calculate_fluxes();
  // TODO: rest

  std::cout << d_Q_up << std::endl;
  std::cout << d_A_up << std::endl;
}

void ExplicitNonlinearFlowSolver::calculate_fluxes() {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q0 = 0;

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
      // TODO: Make this more generic!
      const bool is_inflow_boundary = (vertex->get_coordinate() - lm::Point(0, 0, 0)).norm() < 1e-10;

      const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[0]);

      fe_r.reinit(*edge_r);

      d_dof_map->dof_indices(*edge_r, Q_dof_indices_r, 0);
      d_dof_map->dof_indices(*edge_r, A_dof_indices_r, 1);

      extract_dof(Q_dof_indices_r, d_u_prev, Q_prev_loc_r);
      extract_dof(A_dof_indices_r, d_u_prev, A_prev_loc_r);

      fe_r.evaluate_dof_at_quadrature_points(Q_prev_loc_r,  Q_prev_qp_r);
      fe_r.evaluate_dof_at_quadrature_points(A_prev_loc_r,  A_prev_qp_r);

      // inflow boundary
      if (is_inflow_boundary) {
        // we assert that the edge direction fits our assumptions
        // TODO: Make this assumption more generic!
        assert(edge_r->get_vertex_neighbors()[0] == vertex->get_id());

        const double Q = Q_prev_qp_r[0];
        const double A = A_prev_qp_r[0];

        const double W1 = calculate_W1_value(Q, A, d_G0, d_rho, d_A0);
        const double W2 = 2 * d_inflow_value_function(d_t_now) + calculate_W1_value(Q0, d_A0, d_G0, d_rho, d_A0);

        double Q_up = 0, A_up = 0;
        solve_W12(Q_up, A_up, W1, W2, d_G0, d_rho, d_A0);

        d_Q_up[vertex->get_id()] = Q_up;
        d_A_up[vertex->get_id()] = A_up;
      }
      // free outflow boundary
      else {
        // we assert that the edge direction fits our assumptions
        // TODO: Make this assumption more generic!
        assert(edge_r->get_vertex_neighbors()[1] == vertex->get_id());

        const double Q = Q_prev_qp_r[1];
        const double A = A_prev_qp_r[1];

        const double W1 = calculate_W1_value(Q0, d_A0, d_G0, d_rho, d_A0);
        const double W2 = calculate_W2_value(Q, A, d_G0, d_rho, d_A0);

        double Q_up = 0, A_up = 0;
        solve_W12(Q_up, A_up, W1, W2, d_G0, d_rho, d_A0);

        d_Q_up[vertex->get_id()] = Q_up;
        d_A_up[vertex->get_id()] = A_up;
      }
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

void ExplicitNonlinearFlowSolver::calculate_rhs() {
  // first zero rhs
  for (std::size_t idx=0; idx<d_rhs.size(); idx+=1)
    d_rhs[idx] = 0;

  const std::size_t num_basis_functions = degree + 1;

  // assemble cell contributions
  {
    std::vector<double> f_loc_Q(num_basis_functions);
    std::vector<double> f_loc_A(num_basis_functions);

    FETypeNetwork<degree> fe(create_gauss4());
    const auto &phi = fe.get_phi();
    const auto &dphi = fe.get_dphi();
    const auto &JxW = fe.get_JxW();

    std::vector<std::size_t> Q_dof_indices;
    std::vector<std::size_t> A_dof_indices;

    std::vector<double> Q_prev_loc;
    std::vector<double> A_prev_loc;

    std::vector<double> Q_prev_qp;
    std::vector<double> A_prev_qp;

    for (const auto &e_id : d_graph->get_edge_ids()) {
      const auto edge = d_graph->get_edge(e_id);
      fe.reinit(*edge);

      d_dof_map->dof_indices(*edge, Q_dof_indices, 0);
      d_dof_map->dof_indices(*edge, A_dof_indices, 1);

      extract_dof(Q_dof_indices, d_u_prev, Q_prev_loc);
      extract_dof(A_dof_indices, d_u_prev, A_prev_loc);

      fe.evaluate_dof_at_quadrature_points(Q_prev_loc,  Q_prev_qp);
      fe.evaluate_dof_at_quadrature_points(A_prev_loc,  A_prev_qp);

      // evaluate F = (F_Q, F_A) at the quadrature points
      const auto& F_Q = Q_prev_qp;
      std::vector<double> F_A (Q_prev_qp.size(), 0);
      for (std::size_t qp = 0; qp < phi[0].size(); qp += 1)
        F_A[qp] = std::pow(Q_prev_qp[qp], 2)/A_prev_qp[qp] + d_G0/(3*d_rho*std::sqrt(d_A0)) * std::pow(A_prev_qp[qp], 3./2.);

      // evaluate S = (0, S_A) at the quadrature points
      const double S_Q = 0;
      std::vector<double> S_A (Q_prev_qp.size(), 0);
      for (std::size_t qp = 0; qp < phi[0].size(); qp += 1)
        S_A[qp] = - 2 * d_mu * M_PI * (d_gamma + 2) * Q_prev_qp[qp] / A_prev_qp[qp];

      for (std::size_t i = 0; i < num_basis_functions; i += 1) {
        // rhs integral
        f_loc_A[i] = 0;
        f_loc_Q[i] = 0;
        for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
          // mass term
          f_loc_A[i] += phi[i][qp] * A_prev_qp[qp] * JxW[qp];
          f_loc_A[i] += phi[i][qp] * d_tau * S_A[qp] * JxW[qp];
          f_loc_A[i] += dphi[i][qp] * d_tau * F_A[qp] * JxW[qp];

          f_loc_Q[i] += phi[i][qp] * Q_prev_qp[qp] * JxW[qp];
          f_loc_Q[i] += phi[i][qp] * d_tau * S_Q * JxW[qp];
          f_loc_Q[i] += dphi[i][qp] * d_tau * F_Q[qp] * JxW[qp];
        }
      }

      // copy into global vector
      for (std::size_t i = 0; i < Q_dof_indices.size(); i += 1)
        d_rhs[Q_dof_indices[i]] += f_loc_Q[i];
      for (std::size_t i = 0; i < A_dof_indices.size(); i += 1)
        d_rhs[A_dof_indices[i]] += f_loc_A[i];
    }
  }

  // assemble boundary contributions
  {
    // functions evaluated on the inner cell boundaries
    FETypeInnerBdryNetwork<degree> fe_inner;
    const auto &phi_l = fe_inner.get_phi_l();
    const auto &phi_r = fe_inner.get_phi_r();

    // block rhs for inner boundaries
    std::vector<double> f_loc_l_Q(num_basis_functions);
    std::vector<double> f_loc_l_A(num_basis_functions);
    std::vector<double> f_loc_r_Q(num_basis_functions);
    std::vector<double> f_loc_r_A(num_basis_functions);

    // functions evaluated on the exterior boundaries
    FETypeExteriorBdryNetwork<degree> fe_ext;
    const auto &phi = fe_ext.get_phi();

    // right hand side for inner boundaries
    std::vector<double> f_ext_loc_Q(num_basis_functions);
    std::vector<double> f_ext_loc_A(num_basis_functions);

    std::vector<std::size_t> Q_dof_indices;
    std::vector<std::size_t> A_dof_indices;
    std::vector<std::size_t> Q_dof_indices_l;
    std::vector<std::size_t> A_dof_indices_l;
    std::vector<std::size_t> Q_dof_indices_r;
    std::vector<std::size_t> A_dof_indices_r;

    for (const auto &v_id : d_graph->get_vertex_ids()) {
      const auto vertex = d_graph->get_vertex(v_id);

      // exterior boundary
      if (vertex->is_leaf()) {
        const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
        d_dof_map->dof_indices(*edge, Q_dof_indices, 0);
        d_dof_map->dof_indices(*edge, A_dof_indices, 1);
        fe_ext.reinit(*vertex, *edge);

        // zero local system
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          f_ext_loc_Q[i] = 0;
          f_ext_loc_A[i] = 0;
        }

        // TODO: Make this more generic!
        const bool is_inflow_boundary = (vertex->get_coordinate() - lm::Point(0, 0, 0)).norm() < 1e-10;

        // inflow boundary
        if (is_inflow_boundary) {
          // TODO: implement
        }
        // outflow boundary
        else {
          // TODO: implement
        }

        // copy into vector
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          d_rhs[Q_dof_indices[i]] += f_ext_loc_Q[i];
          d_rhs[A_dof_indices[i]] += f_ext_loc_A[i];
        }
      }
        // inner boundary
      else {
        const auto edge_l = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
        const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[1]);

        d_dof_map->dof_indices(*edge_l, Q_dof_indices_l, 0);
        d_dof_map->dof_indices(*edge_l, A_dof_indices_l, 1);
        d_dof_map->dof_indices(*edge_r, Q_dof_indices_r, 0);
        d_dof_map->dof_indices(*edge_r, A_dof_indices_r, 1);

        fe_inner.reinit(*vertex, *edge_l, *edge_r);

        // zero local system
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          f_loc_l_Q[i] = 0;
          f_loc_l_A[i] = 0;
          f_loc_r_Q[i] = 0;
          f_loc_r_A[i] = 0;
        }

        // TODO: assemble contributions

        // copy into global matrix
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          d_rhs[Q_dof_indices_l[i]] += f_loc_l_Q[i];
          d_rhs[A_dof_indices_l[i]] += f_loc_l_A[i];
          d_rhs[Q_dof_indices_r[i]] += f_loc_r_Q[i];
          d_rhs[A_dof_indices_r[i]] += f_loc_r_A[i];
        }
      }
    }
  }
}

void ExplicitNonlinearFlowSolver::apply_inverse_mass() {
  // TODO: implement
}

} // namespace macrocirculation