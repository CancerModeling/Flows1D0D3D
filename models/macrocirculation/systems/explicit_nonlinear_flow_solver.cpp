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
#include "libmesh/libmesh.h"

namespace macrocirculation {

namespace lm = libMesh;

constexpr std::size_t degree = 1;

void interpolate_constant(const GraphStorage & graph, const DofMapNetwork & dof_map, double value, std::size_t component, std::vector< double > & result)
{
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    dof_map.dof_indices(*edge, dof_indices, component);
    result[dof_indices[0]] = value;
  }
}

void assemble_inverse_mass(const GraphStorage &graph, const DofMapNetwork &dof_map, std::vector<double> &inv_mass) {
  // make sure that the inverse mass vector is large enough
  assert(inv_mass.size() == dof_map.num_dof());

  const std::size_t num_basis_functions = degree + 1;

  // the squared norm of the legendre polynomials on [-1, +1]
  std::vector<double> legendre_weight = {2, 2. / 3., 0.4, 0.285714};

  // assert that we have precalculated enough weights
  assert(num_basis_functions <= legendre_weight.size());

  std::vector<std::size_t> Q_dof_indices(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices(num_basis_functions, 0);

  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);

    dof_map.dof_indices(*edge, Q_dof_indices, 0);
    dof_map.dof_indices(*edge, A_dof_indices, 1);

    const double edge_weight = edge->get_length() / 2;

    for (std::size_t i = 0; i < num_basis_functions; i += 1) {
      const double mass = legendre_weight[i] * edge_weight;
      inv_mass[Q_dof_indices[i]] = 1. / mass;
      inv_mass[A_dof_indices[i]] = 1. / mass;
    }
  }
}

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
      d_inverse_mass(d_dof_map->num_dof()),
      d_G0(592.4e2), // 592.4 10^2 Pa,  TODO: Check if units are consistent!
      d_rho(1.028),  // 1.028 kg/cm^3,  TODO: Check if units are consistent!
      d_A0(6.97),    // 6.97 cm^2,      TODO: Check if units are consistent!
      d_mu(4.5),     // 4.5 m Pa/s      TODO: Check if units are consistent!
      d_gamma(2)     // Poiseuille flow
{
  // set A constant to A0
  interpolate_constant(*d_graph, *d_dof_map, d_A0, 1, d_u_prev);
  interpolate_constant(*d_graph, *d_dof_map, d_A0, 1, d_u_now);
  assemble_inverse_mass(*d_graph, *d_dof_map, d_inverse_mass);
}

void ExplicitNonlinearFlowSolver::solve() {
  calculate_fluxes();
  calculate_rhs();
  apply_inverse_mass();

  std::cout << d_Q_up << std::endl;
  std::cout << d_A_up << std::endl;

  std::cout << d_rhs << std::endl;
  std::cout << d_u_prev << std::endl;
  std::cout << d_u_now << std::endl;
}

void ExplicitNonlinearFlowSolver::calculate_fluxes() {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q_init = 0;
  const double A_init = d_A0;

  const std::size_t num_basis_functions = degree + 1;

  // finite-element for left and right edge
  FETypeNetwork<degree> fe_l(create_trapezoidal_rule());
  FETypeNetwork<degree> fe_r(create_trapezoidal_rule());

  // dof indices for left and right edge
  std::vector<std::size_t> Q_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> Q_dof_indices_r(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_r(num_basis_functions, 0);

  // local views of our previous solution
  std::vector<double> Q_prev_loc_r(num_basis_functions, 0);
  std::vector<double> A_prev_loc_r(num_basis_functions, 0);
  std::vector<double> Q_prev_loc_l(num_basis_functions, 0);
  std::vector<double> A_prev_loc_l(num_basis_functions, 0);

  // previous solution evaluated at quadrature points
  // we have 2 quadrature points for the trapezoidal rule
  std::vector<double> Q_prev_qp_r(2, 0);
  std::vector<double> A_prev_qp_r(2, 0);
  std::vector<double> Q_prev_qp_l(2, 0);
  std::vector<double> A_prev_qp_l(2, 0);

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

      fe_r.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
      fe_r.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

      // inflow boundary
      if (is_inflow_boundary) {
        // we assert that the edge direction fits our assumptions
        // TODO: Make this assumption more generic!
        assert(edge_r->get_vertex_neighbors()[0] == vertex->get_id());

        const double Q = Q_prev_qp_r[0];
        const double A = A_prev_qp_r[0];

        const double W1 = calculate_W1_value(Q, A, d_G0, d_rho, d_A0);
        const double W2 = 2 * d_inflow_value_function(d_t_now) + calculate_W1_value(Q_init, A_init, d_G0, d_rho, d_A0);

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

        const double W1 = calculate_W1_value(Q_init, A_init, d_G0, d_rho, d_A0);
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

      fe_l.evaluate_dof_at_quadrature_points(Q_prev_loc_l, Q_prev_qp_l);
      fe_l.evaluate_dof_at_quadrature_points(A_prev_loc_l, A_prev_qp_l);
      fe_l.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
      fe_l.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

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
  for (std::size_t idx = 0; idx < d_rhs.size(); idx += 1)
    d_rhs[idx] = 0;

  const std::size_t num_basis_functions = degree + 1;

  std::vector<double> f_loc_Q(num_basis_functions);
  std::vector<double> f_loc_A(num_basis_functions);

  // data structures for cell contributions
  FETypeNetwork<degree> fe(create_gauss4());
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  std::vector<std::size_t> Q_dof_indices(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices(num_basis_functions, 0);

  std::vector<double> Q_prev_loc(num_basis_functions, 0);
  std::vector<double> A_prev_loc(num_basis_functions, 0);

  std::vector<double> Q_prev_qp(fe.num_quad_points(), 0);
  std::vector<double> A_prev_qp(fe.num_quad_points(), 0);

  // data structures for facet contributions
  FETypeNetwork<degree> fe_boundary(create_trapezoidal_rule());
  const auto& phi_b = fe_boundary.get_phi();

  for (const auto &e_id : d_graph->get_edge_ids()) {
    const auto edge = d_graph->get_edge(e_id);
    fe.reinit(*edge);
    fe_boundary.reinit(*edge);

    // evaluate Q and A inside cell
    d_dof_map->dof_indices(*edge, Q_dof_indices, 0);
    d_dof_map->dof_indices(*edge, A_dof_indices, 1);

    extract_dof(Q_dof_indices, d_u_prev, Q_prev_loc);
    extract_dof(A_dof_indices, d_u_prev, A_prev_loc);

    fe.evaluate_dof_at_quadrature_points(Q_prev_loc, Q_prev_qp);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc, A_prev_qp);

    // evaluate Q and A on boundary
    const double Q_up_0 = d_Q_up[edge->get_vertex_neighbors()[0]];
    const double Q_up_1 = d_Q_up[edge->get_vertex_neighbors()[1]];
    const double A_up_0 = d_A_up[edge->get_vertex_neighbors()[0]];
    const double A_up_1 = d_A_up[edge->get_vertex_neighbors()[1]];

    // the A-component of our F function
    const auto F_A_eval = [=](double Q, double A) -> double {
      return std::pow(Q, 2) / A + d_G0 / (3 * d_rho * std::sqrt(d_A0)) * std::pow(A, 3. / 2.);
    };

    // evaluate F = (F_Q, F_A) at the quadrature points
    const auto &F_Q = Q_prev_qp;
    std::vector<double> F_A(Q_prev_qp.size(), 0);
    for (std::size_t qp = 0; qp < fe.num_quad_points(); qp += 1)
      F_A[qp] = F_A_eval(Q_prev_qp[qp], A_prev_qp[qp]);

    // evaluate S = (0, S_A) at the quadrature points
    const double S_Q = 0;
    std::vector<double> S_A(Q_prev_qp.size(), 0);
    for (std::size_t qp = 0; qp < fe.num_quad_points(); qp += 1)
      S_A[qp] = -2 * d_mu * M_PI * (d_gamma + 2) * Q_prev_qp[qp] / A_prev_qp[qp];

    for (std::size_t i = 0; i < num_basis_functions; i += 1) {
      // rhs integral
      f_loc_A[i] = 0;
      f_loc_Q[i] = 0;

      // cell contributions
      for (std::size_t qp = 0; qp < fe.num_quad_points(); qp += 1) {
        f_loc_Q[i] += phi[i][qp] * Q_prev_qp[qp] * JxW[qp];
        f_loc_Q[i] += phi[i][qp] * d_tau * S_Q * JxW[qp];
        f_loc_Q[i] += dphi[i][qp] * d_tau * F_Q[qp] * JxW[qp];

        f_loc_A[i] += phi[i][qp] * A_prev_qp[qp] * JxW[qp];
        f_loc_A[i] += phi[i][qp] * d_tau * S_A[qp] * JxW[qp];
        f_loc_A[i] += dphi[i][qp] * d_tau * F_A[qp] * JxW[qp];
      }

      // boundary contributions
      f_loc_Q[i] += d_tau * Q_up_1 * phi_b[i][1];
      f_loc_Q[i] -= d_tau * Q_up_0 * phi_b[i][0];

      f_loc_A[i] += d_tau * F_A_eval(Q_up_1, A_up_1) * phi_b[i][1];
      f_loc_A[i] -= d_tau * F_A_eval(Q_up_0, A_up_0) * phi_b[i][0];
    }

    // copy into global vector
    for (std::size_t i = 0; i < Q_dof_indices.size(); i += 1)
      d_rhs[Q_dof_indices[i]] += f_loc_Q[i];
    for (std::size_t i = 0; i < A_dof_indices.size(); i += 1)
      d_rhs[A_dof_indices[i]] += f_loc_A[i];
  }
}

void ExplicitNonlinearFlowSolver::apply_inverse_mass() {
  for (std::size_t i=0; i<d_dof_map->num_dof(); i+=1)
    d_u_now[i] = d_inverse_mass[i] * d_rhs[i];
}

} // namespace macrocirculation