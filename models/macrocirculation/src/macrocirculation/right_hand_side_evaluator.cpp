////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "right_hand_side_evaluator.hpp"

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "libmesh/point.h"
#include "vessel_data_storage.hpp"
#include "vessel_formulas.hpp"
#include <utility>
#include "communication/mpi.hpp"
#include "communicator.hpp"

constexpr const std::size_t degree = 1;

namespace macrocirculation {

namespace lm = libMesh;

default_S::default_S(double mu, double gamma, double phi)
    : d_mu(mu), d_gamma(gamma), d_phi(phi) {}

void default_S::operator()(double, const std::vector<Point>&, const std::vector<double> &Q, const std::vector<double> &A, std::vector<double> &S_Q_out, std::vector<double> &S_A_out) const {
  // all vectors have to have the same shape
  assert(Q.size() == A.size());
  assert(Q.size() == S_Q_out.size());
  assert(Q.size() == S_A_out.size());

  for (std::size_t qp = 0; qp < Q.size(); qp += 1) {
    S_Q_out[qp] = -2 * d_mu * M_PI * (d_gamma + 2) * Q[qp] / A[qp];
    S_A_out[qp] = d_phi;
  }
}

template<std::size_t degree>
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

template<std::size_t degree>
RightHandSideEvaluator<degree>::RightHandSideEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data, std::shared_ptr<DofMapNetwork> dof_map)
    : d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_vessel_data(std::move(vessel_data)),
      d_communicator(std::make_unique<Communicator>(comm, d_graph, d_dof_map)),
      d_S_evaluator(default_S{
        4.5,    // 4.5 m Pa/s
        2,   // Poiseuille flow
        0       // 0 cm^2/s, no wall permeability
      }),
      d_Q_up_el(d_graph->num_edges()),
      d_Q_up_er(d_graph->num_edges()),
      d_A_up_el(d_graph->num_edges()),
      d_A_up_er(d_graph->num_edges()),
      d_inverse_mass(d_dof_map->num_dof()),
      d_comm(comm) {
  assemble_inverse_mass<degree>(*d_graph, *d_dof_map, d_inverse_mass);
}

template<std::size_t degree>
void RightHandSideEvaluator<degree>::evaluate(const double t, std::vector<double> &u_prev, std::vector<double> &rhs) {
  //
  d_communicator->update_ghost_layer(u_prev);
  calculate_fluxes(t, u_prev);
  calculate_rhs(t, u_prev, rhs);
  apply_inverse_mass(rhs);
}

template<std::size_t degree>
void RightHandSideEvaluator<degree>::set_rhs_S(VectorEvaluator S_evaluator) { d_S_evaluator = S_evaluator; }

template<std::size_t degree>
void RightHandSideEvaluator<degree>::calculate_fluxes(const double t, const std::vector<double> &u_prev) {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q_init = 0;

  const std::size_t num_basis_functions = degree + 1;

  // finite-element for left and right edge
  FETypeNetwork<degree> fe_l(create_trapezoidal_rule());
  FETypeNetwork<degree> fe_r(create_trapezoidal_rule());
  FETypeExteriorBdryNetwork<degree> fe_ext;

  // dof indices for left and right edge
  std::vector<std::size_t> Q_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> Q_dof_indices_r(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_r(num_basis_functions, 0);

  // local views of our previous solution
  std::vector<double> Q_prev_loc_l(num_basis_functions, 0);
  std::vector<double> A_prev_loc_l(num_basis_functions, 0);
  std::vector<double> Q_prev_loc_r(num_basis_functions, 0);
  std::vector<double> A_prev_loc_r(num_basis_functions, 0);

  // previous solution evaluated at quadrature points
  // we have 2 quadrature points for the trapezoidal rule
  std::vector<double> Q_prev_qp_l(2, 0);
  std::vector<double> A_prev_qp_l(2, 0);
  std::vector<double> Q_prev_qp_r(2, 0);
  std::vector<double> A_prev_qp_r(2, 0);

  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // exterior boundary
    if (vertex->is_leaf()) {
      const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[0]);

      const auto &param = d_vessel_data->get_parameters(*edge_r);

      fe_r.reinit(*edge_r);

      d_dof_map->dof_indices(*edge_r, Q_dof_indices_r, 0);
      d_dof_map->dof_indices(*edge_r, A_dof_indices_r, 1);

      extract_dof(Q_dof_indices_r, u_prev, Q_prev_loc_r);
      extract_dof(A_dof_indices_r, u_prev, A_prev_loc_r);

      fe_r.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
      fe_r.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

      // inflow boundary
      if (vertex->is_inflow()) {
        // does the vessel point towards the vertex?
        const bool in = edge_r->get_vertex_neighbors()[1] == vertex->get_id();

        const double Q = in ? Q_prev_qp_r[1] : Q_prev_qp_r[0];
        const double A = in ? A_prev_qp_r[1] : A_prev_qp_r[0];

        const double Q_star = vertex->get_inflow_value(t);
        const double A_up = assemble_in_flow(Q, A, in, Q_star, param.G0, param.rho, param.A0);

        if (in)
        {
          d_Q_up_er[edge_r->get_id()] = Q_star;
          d_A_up_er[edge_r->get_id()] = A_up;
        }
        else
        {
          d_Q_up_el[edge_r->get_id()] = Q_star;
          d_A_up_el[edge_r->get_id()] = A_up;
        }
      }
      // free outflow boundary
      else {
        // we assert that the edge direction fits our assumptions
        // TODO: Make this assumption more generic!
        assert(edge_r->get_vertex_neighbors()[1] == vertex->get_id());

        // TODO: make this more generic for other initial flow values
        const double A_init = param.A0;

        const double Q = Q_prev_qp_r[1];
        const double A = A_prev_qp_r[1];

        const double W1 = calculate_W1_value(Q_init, A_init, param.G0, param.rho, param.A0);
        const double W2 = calculate_W2_value(Q, A, param.G0, param.rho, param.A0);

        double Q_up = 0, A_up = 0;
        solve_W12(Q_up, A_up, W1, W2, param.G0, param.rho, param.A0);

        d_Q_up_er[edge_r->get_id()] = Q_up;
        d_A_up_er[edge_r->get_id()] = A_up;
      }
    }
    // bifurcation boundary
    else if (vertex->is_bifurcation()) {
      if (vertex->get_edge_neighbors().size() > 3)
        throw std::runtime_error("only 3 neighbors at bifurcation points possible.");

      // get vertices
      const auto e0 = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      const auto e1 = d_graph->get_edge(vertex->get_edge_neighbors()[1]);
      const auto e2 = d_graph->get_edge(vertex->get_edge_neighbors()[2]);

      // check orientation
      const bool e0_in = e0->is_pointing_to(vertex->get_id());
      const bool e1_in = e1->is_pointing_to(vertex->get_id());
      const bool e2_in = e2->is_pointing_to(vertex->get_id());

      // get data
      const auto &p_e0 = d_vessel_data->get_parameters(*e0);
      const auto &p_e1 = d_vessel_data->get_parameters(*e1);
      const auto &p_e2 = d_vessel_data->get_parameters(*e2);

      // evaluate on edge e0
      d_dof_map->dof_indices(*e0, Q_dof_indices_l, 0);
      d_dof_map->dof_indices(*e0, A_dof_indices_l, 1);
      extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
      extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
      fe_ext.reinit(*vertex, *e0);
      const double Q_e0 = fe_ext.evaluate_dof_at_boundary_points(Q_prev_loc_l);
      const double A_e0 = fe_ext.evaluate_dof_at_boundary_points(A_prev_loc_l);

      // evaluate on edge e1
      d_dof_map->dof_indices(*e1, Q_dof_indices_l, 0);
      d_dof_map->dof_indices(*e1, A_dof_indices_l, 1);
      extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
      extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
      fe_ext.reinit(*vertex, *e1);
      const double Q_e1 = fe_ext.evaluate_dof_at_boundary_points(Q_prev_loc_l);
      const double A_e1 = fe_ext.evaluate_dof_at_boundary_points(A_prev_loc_l);

      // evaluate on edge e2
      d_dof_map->dof_indices(*e2, Q_dof_indices_l, 0);
      d_dof_map->dof_indices(*e2, A_dof_indices_l, 1);
      extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
      extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
      fe_ext.reinit(*vertex, *e2);
      const double Q_e2 = fe_ext.evaluate_dof_at_boundary_points(Q_prev_loc_l);
      const double A_e2 = fe_ext.evaluate_dof_at_boundary_points(A_prev_loc_l);

      // get upwinded values at bifurcation
      double Q_e0_up, A_e0_up;
      double Q_e1_up, A_e1_up;
      double Q_e2_up, A_e2_up;
      const auto num_iter = solve_at_bifurcation(
        Q_e0, A_e0, p_e0, e0_in,
        Q_e1, A_e1, p_e1, e1_in,
        Q_e2, A_e2, p_e2, e2_in,
        Q_e0_up, A_e0_up,
        Q_e1_up, A_e1_up,
        Q_e2_up, A_e2_up);

      // save upwinded values into upwind vector
      if (e0_in) {
        d_Q_up_er[e0->get_id()] = Q_e0_up;
        d_A_up_er[e0->get_id()] = A_e0_up;
      } else {
        d_Q_up_el[e0->get_id()] = Q_e0_up;
        d_A_up_el[e0->get_id()] = A_e0_up;
      }

      if (e1_in) {
        d_Q_up_er[e1->get_id()] = Q_e1_up;
        d_A_up_er[e1->get_id()] = A_e1_up;
      } else {
        d_Q_up_el[e1->get_id()] = Q_e1_up;
        d_A_up_el[e1->get_id()] = A_e1_up;
      }

      if (e2_in) {
        d_Q_up_er[e2->get_id()] = Q_e2_up;
        d_A_up_er[e2->get_id()] = A_e2_up;
      } else {
        d_Q_up_el[e2->get_id()] = Q_e2_up;
        d_A_up_el[e2->get_id()] = A_e2_up;
      }
    }
    // inner boundary
    else {
      const auto edge_l = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[1]);

      const auto &p_l = d_vessel_data->get_parameters(*edge_l);
      const auto &p_r = d_vessel_data->get_parameters(*edge_r);

      // we assert that the edge direction fits our assumptions
      assert(edge_l->get_vertex_neighbors()[1] == vertex->get_id());
      assert(edge_r->get_vertex_neighbors()[0] == vertex->get_id());

      fe_l.reinit(*edge_l);
      fe_r.reinit(*edge_r);

      d_dof_map->dof_indices(*edge_l, Q_dof_indices_l, 0);
      d_dof_map->dof_indices(*edge_l, A_dof_indices_l, 1);
      d_dof_map->dof_indices(*edge_r, Q_dof_indices_r, 0);
      d_dof_map->dof_indices(*edge_r, A_dof_indices_r, 1);

      extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
      extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
      extract_dof(Q_dof_indices_r, u_prev, Q_prev_loc_r);
      extract_dof(A_dof_indices_r, u_prev, A_prev_loc_r);

      fe_l.evaluate_dof_at_quadrature_points(Q_prev_loc_l, Q_prev_qp_l);
      fe_l.evaluate_dof_at_quadrature_points(A_prev_loc_l, A_prev_qp_l);
      fe_r.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
      fe_r.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

      const double Q_l = Q_prev_qp_l[1];
      const double A_l = A_prev_qp_l[1];
      const double Q_r = Q_prev_qp_r[0];
      const double A_r = A_prev_qp_r[0];

      const double W2_l = calculate_W2_value(Q_l, A_l, p_l.G0, p_l.rho, p_l.A0);
      const double W1_r = calculate_W1_value(Q_r, A_r, p_r.G0, p_r.rho, p_r.A0);

      double Q_up = 0, A_up = 0;

      // TODO: make more general, such that different parameters are possible
      assert(p_l.G0 == p_r.G0 && p_l.rho == p_r.rho && p_l.A0 == p_r.A0);
      solve_W12(Q_up, A_up, W1_r, W2_l, p_l.G0, p_l.rho, p_l.A0);

      d_Q_up_er[edge_l->get_id()] = Q_up;
      d_A_up_er[edge_l->get_id()] = A_up;
      d_Q_up_el[edge_r->get_id()] = Q_up;
      d_A_up_el[edge_r->get_id()] = A_up;
    }
  }
}

template<std::size_t degree>
void RightHandSideEvaluator<degree>::calculate_rhs(const double t, const std::vector<double> &u_prev, std::vector<double> &rhs) {
  // first zero rhs
  for (std::size_t idx = 0; idx < rhs.size(); idx += 1)
    rhs[idx] = 0;

  const std::size_t num_basis_functions = degree + 1;

  std::vector<double> f_loc_Q(num_basis_functions);
  std::vector<double> f_loc_A(num_basis_functions);

  // data structures for cell contributions
  QuadratureFormula qf = create_gauss4();
  FETypeNetwork<degree> fe(qf);
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  // some right hand sides need the physical location of the quadrature points
  QuadraturePointMapper qpm(qf);
  const auto &points = qpm.get_quadrature_points();

  std::vector<std::size_t> Q_dof_indices(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices(num_basis_functions, 0);

  std::vector<double> Q_prev_loc(num_basis_functions, 0);
  std::vector<double> A_prev_loc(num_basis_functions, 0);

  std::vector<double> Q_prev_qp(fe.num_quad_points(), 0);
  std::vector<double> A_prev_qp(fe.num_quad_points(), 0);

  const auto &F_A = Q_prev_qp;
  std::vector<double> F_Q(Q_prev_qp.size(), 0);

  std::vector<double> S_Q(fe.num_quad_points(), 0);
  std::vector<double> S_A(fe.num_quad_points(), 0);

  // data structures for facet contributions
  FETypeNetwork<degree> fe_boundary(create_trapezoidal_rule());
  const auto &phi_b = fe_boundary.get_phi();

  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto edge = d_graph->get_edge(e_id);
    fe.reinit(*edge);
    fe_boundary.reinit(*edge);
    qpm.reinit(*edge);

    const auto &param = d_vessel_data->get_parameters(*edge);

    // evaluate Q and A inside cell
    d_dof_map->dof_indices(*edge, Q_dof_indices, 0);
    d_dof_map->dof_indices(*edge, A_dof_indices, 1);

    extract_dof(Q_dof_indices, u_prev, Q_prev_loc);
    extract_dof(A_dof_indices, u_prev, A_prev_loc);

    fe.evaluate_dof_at_quadrature_points(Q_prev_loc, Q_prev_qp);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc, A_prev_qp);

    // evaluate Q and A on boundary
    const double Q_up_0 = d_Q_up_el[edge->get_id()];
    const double Q_up_1 = d_Q_up_er[edge->get_id()];
    const double A_up_0 = d_A_up_el[edge->get_id()];
    const double A_up_1 = d_A_up_er[edge->get_id()];

    // the A-component of our F function
    const auto F_Q_eval = [&param](double Q, double A) -> double {
      return std::pow(Q, 2) / A + param.G0 / (3 * param.rho * std::sqrt(param.A0)) * std::pow(A, 3. / 2.);
    };

    // evaluate F = (F_Q, F_A) at the quadrature points
    for (std::size_t qp = 0; qp < fe.num_quad_points(); qp += 1)
      F_Q[qp] = F_Q_eval(Q_prev_qp[qp], A_prev_qp[qp]);

    // evaluate S = (S_Q, S_A) at the quadrature points
    d_S_evaluator(t, points, Q_prev_qp, A_prev_qp, S_Q, S_A);

    for (std::size_t i = 0; i < num_basis_functions; i += 1) {
      // rhs integral
      f_loc_A[i] = 0;
      f_loc_Q[i] = 0;

      // cell contributions
      for (std::size_t qp = 0; qp < fe.num_quad_points(); qp += 1) {
        f_loc_Q[i] += phi[i][qp] * S_Q[qp] * JxW[qp];
        f_loc_Q[i] += dphi[i][qp] * F_Q[qp] * JxW[qp];

        f_loc_A[i] += phi[i][qp] * S_A[qp] * JxW[qp];
        f_loc_A[i] += dphi[i][qp] * F_A[qp] * JxW[qp];
      }

      // boundary contributions  - tau [ F(U_up) phi ] ds, keep attention to the minus!
      f_loc_Q[i] -= F_Q_eval(Q_up_1, A_up_1) * phi_b[i][1];
      f_loc_Q[i] += F_Q_eval(Q_up_0, A_up_0) * phi_b[i][0];

      f_loc_A[i] -= Q_up_1 * phi_b[i][1];
      f_loc_A[i] += Q_up_0 * phi_b[i][0];
    }

    // copy into global vector
    for (std::size_t i = 0; i < Q_dof_indices.size(); i += 1)
      rhs[Q_dof_indices[i]] += f_loc_Q[i];
    for (std::size_t i = 0; i < A_dof_indices.size(); i += 1)
      rhs[A_dof_indices[i]] += f_loc_A[i];
  }
}

template<std::size_t degree>
void RightHandSideEvaluator<degree>::apply_inverse_mass(std::vector<double> &rhs) {
  for (std::size_t i = 0; i < d_dof_map->num_dof(); i += 1)
    rhs[i] = d_inverse_mass[i] * rhs[i];
}

template void assemble_inverse_mass<0>(const GraphStorage &, const DofMapNetwork &, std::vector<double> &);
template void assemble_inverse_mass<1>(const GraphStorage &, const DofMapNetwork &, std::vector<double> &);
template void assemble_inverse_mass<2>(const GraphStorage &, const DofMapNetwork &, std::vector<double> &);
template void assemble_inverse_mass<3>(const GraphStorage &, const DofMapNetwork &, std::vector<double> &);

template class RightHandSideEvaluator<0>;
template class RightHandSideEvaluator<1>;
template class RightHandSideEvaluator<2>;
template class RightHandSideEvaluator<3>;

} // namespace macrocirculation