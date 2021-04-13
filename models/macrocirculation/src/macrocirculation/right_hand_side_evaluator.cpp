////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "right_hand_side_evaluator.hpp"

#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

default_S::default_S(double mu, double gamma, double phi)
    : d_mu(mu),
      d_gamma(gamma),
      d_phi(phi) {}

void default_S::operator()(double,
                           const std::vector<double> &,
                           const std::vector<double> &Q,
                           const std::vector<double> &A,
                           std::vector<double> &S_Q_out,
                           std::vector<double> &S_A_out) const {
  // all vectors have to have the same shape
  assert(Q.size() == A.size());
  assert(Q.size() == S_Q_out.size());
  assert(Q.size() == S_A_out.size());

  for (std::size_t qp = 0; qp < Q.size(); qp += 1) {
    S_Q_out[qp] = -2 * d_mu * M_PI * (d_gamma + 2) * Q[qp] / A[qp];
    S_A_out[qp] = d_phi;
  }
}

void assemble_inverse_mass(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, std::vector<double> &inv_mass) {
  // make sure that the inverse mass vector is large enough
  assert(inv_mass.size() == dof_map.num_dof());

  // the squared norm of the legendre polynomials on [-1, +1]
  std::vector<double> legendre_weight = {2, 2. / 3., 0.4, 0.285714};

  std::vector<std::size_t> dof_indices(4, 0);

  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);
    const auto local_dof_map = dof_map.get_local_dof_map(*edge);

    const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

    // assert that we have precalculated enough weights
    assert(num_basis_functions <= legendre_weight.size());

    dof_indices.resize(num_basis_functions);

    const auto edge_length = edge->get_physical_data().length;
    const auto micro_edge_length = edge_length / local_dof_map.num_micro_edges();
    const double edge_weight = micro_edge_length / 2;

    for (std::size_t component = 0; component < local_dof_map.num_components(); component += 1) {
      for (std::size_t local_micro_edge_id = 0; local_micro_edge_id < local_dof_map.num_micro_edges();
           local_micro_edge_id += 1) {
        local_dof_map.dof_indices(local_micro_edge_id, component, dof_indices);

        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          const double mass = legendre_weight[i] * edge_weight;
          inv_mass[dof_indices[i]] = 1. / mass;
        }
      }
    }
  }
}

RightHandSideEvaluator::RightHandSideEvaluator(MPI_Comm comm,
                                               std::shared_ptr<GraphStorage> graph,
                                               std::shared_ptr<DofMap> dof_map,
                                               std::size_t degree)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_edge_boundary_communicator(Communicator::create_edge_boundary_value_communicator(comm, d_graph)),
      d_S_evaluator(default_S{
        4.5, // 4.5 m Pa/s
        9,   // Poiseuille flow
        0    // 0 cm^2/s, no wall permeability
      }),
      d_degree(degree),
      d_Q_macro_edge_boundary_value(2 * d_graph->num_edges()),
      d_A_macro_edge_boundary_value(2 * d_graph->num_edges()),
      d_Q_macro_edge_flux_l(d_graph->num_edges()),
      d_Q_macro_edge_flux_r(d_graph->num_edges()),
      d_A_macro_edge_flux_l(d_graph->num_edges()),
      d_A_macro_edge_flux_r(d_graph->num_edges()),
      d_inverse_mass(d_dof_map->num_dof()) {
  assemble_inverse_mass(d_comm, *d_graph, *d_dof_map, d_inverse_mass);
}

void RightHandSideEvaluator::evaluate(const double t, const std::vector<double> &u_prev, std::vector<double> &rhs) {
  evaluate_macro_edge_boundary_values(u_prev);
  d_edge_boundary_communicator.update_ghost_layer(d_Q_macro_edge_boundary_value);
  d_edge_boundary_communicator.update_ghost_layer(d_A_macro_edge_boundary_value);
  //std::cout << "fluxes Q_macro_edge_boundary_value_l " << d_Q_macro_edge_boundary_value_l << std::endl;
  //std::cout << "fluxes Q_macro_edge_boundary_value_r " << d_Q_macro_edge_boundary_value_r << std::endl;
  //std::cout << "fluxes A_macro_edge_boundary_value_l " << d_A_macro_edge_boundary_value_l << std::endl;
  //std::cout << "fluxes A_macro_edge_boundary_value_r " << d_A_macro_edge_boundary_value_r << std::endl;
  calculate_fluxes(t, u_prev);
  //std::cout << "fluxes d_Q_macro_edge_flux_r " << d_Q_macro_edge_flux_r << std::endl;
  //std::cout << "fluxes d_A_macro_edge_flux_r " << d_A_macro_edge_flux_r << std::endl;
  //std::cout << "fluxes d_Q_macro_edge_flux_l " << d_Q_macro_edge_flux_l << std::endl;
  //std::cout << "fluxes d_A_macro_edge_flux_l " << d_A_macro_edge_flux_l << std::endl;
  calculate_rhs(t, u_prev, rhs);
  //std::cout << "rhs " << rhs << std::endl;
  apply_inverse_mass(rhs);
  if (std::isnan(rhs.front()) || std::isnan(rhs.back())) {
    throw std::runtime_error("contains nans");
  }
}

void RightHandSideEvaluator::set_rhs_S(VectorEvaluator S_evaluator) {
  d_S_evaluator = std::move(S_evaluator);
}

void RightHandSideEvaluator::evaluate_macro_edge_boundary_values(const std::vector<double> &u_prev) {
  std::vector<std::size_t> dof_indices(4, 0);
  std::vector<double> local_dofs(4, 0);

  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto edge = d_graph->get_edge(e_id);
    const auto &param = edge->get_physical_data();
    const auto &local_dof_map = d_dof_map->get_local_dof_map(*edge);
    const double h = param.length / local_dof_map.num_micro_edges();

    FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    dof_indices.resize(local_dof_map.num_basis_functions());
    local_dofs.resize(local_dof_map.num_basis_functions());

    local_dof_map.dof_indices(0, 0, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_Q_macro_edge_boundary_value[2 * edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

    local_dof_map.dof_indices(0, 1, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_A_macro_edge_boundary_value[2 * edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

    local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 0, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_Q_macro_edge_boundary_value[2 * edge->get_id() + 1] = fe.evaluate_dof_at_boundary_points(local_dofs).right;

    local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 1, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_A_macro_edge_boundary_value[2 * edge->get_id() + 1] = fe.evaluate_dof_at_boundary_points(local_dofs).right;
  }
}

void RightHandSideEvaluator::calculate_nfurcation_fluxes(const std::vector<double> &u_prev) {
  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // we only handle bifurcations
    if (!vertex->is_bifurcation())
      continue;

    const size_t num_vessels = vertex->get_edge_neighbors().size();

    // get edges
    std::vector<std::shared_ptr<Edge>> e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      e.push_back(d_graph->get_edge(vertex->get_edge_neighbors()[vessel_idx]));

    // check orientation
    std::vector<bool> e_in;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      e_in.push_back(e[vessel_idx]->is_pointing_to(vertex->get_id()));

    // get data
    std::vector<VesselParameters> p_e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      const auto &data_e = e[vessel_idx]->get_physical_data();
      p_e.emplace_back(data_e.G0, data_e.A0, data_e.rho);
    }

    // evaluate on edges
    std::vector<double> Q_e;
    std::vector<double> A_e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      const double Q = e_in[vessel_idx] ? d_Q_macro_edge_boundary_value[2 * e[vessel_idx]->get_id() + 1] : d_Q_macro_edge_boundary_value[2 * e[vessel_idx]->get_id()];
      const double A = e_in[vessel_idx] ? d_A_macro_edge_boundary_value[2 * e[vessel_idx]->get_id() + 1] : d_A_macro_edge_boundary_value[2 * e[vessel_idx]->get_id()];
      Q_e.push_back(Q);
      A_e.push_back(A);
    }

    // get upwinded values at bifurcation
    std::vector<double> Q_up(num_vessels, 0);
    std::vector<double> A_up(num_vessels, 0);
    solve_at_nfurcation(Q_e, A_e, p_e, e_in, Q_up, A_up);

    // save upwinded values into upwind vector
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      if (e_in[vessel_idx]) {
        d_Q_macro_edge_flux_r[e[vessel_idx]->get_id()] = Q_up[vessel_idx];
        d_A_macro_edge_flux_r[e[vessel_idx]->get_id()] = A_up[vessel_idx];
      } else {
        d_Q_macro_edge_flux_l[e[vessel_idx]->get_id()] = Q_up[vessel_idx];
        d_A_macro_edge_flux_l[e[vessel_idx]->get_id()] = A_up[vessel_idx];
      }
    }
  }
}

void RightHandSideEvaluator::calculate_inout_fluxes(double t, const std::vector<double> &u_prev) {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q_init = 0;

  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // exterior boundary
    if (vertex->is_leaf()) {
      const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);

      const auto &param = edge->get_physical_data();

      // does the vessel point towards the vertex?
      const bool in = edge->is_pointing_to(vertex->get_id());

      const double Q =
        in ? d_Q_macro_edge_boundary_value[2 * edge->get_id() + 1] : d_Q_macro_edge_boundary_value[2 * edge->get_id()];
      const double A =
        in ? d_A_macro_edge_boundary_value[2 * edge->get_id() + 1] : d_A_macro_edge_boundary_value[2 * edge->get_id()];

      // inflow boundary
      if (vertex->is_inflow()) {
        const double Q_star = (in ? -1 : + 1) * vertex->get_inflow_value(t);
        const double A_up = assemble_in_flow(Q, A, in, Q_star, param.G0, param.rho, param.A0);

        if (in) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_star;
          d_A_macro_edge_flux_r[edge->get_id()] = A_up;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_star;
          d_A_macro_edge_flux_l[edge->get_id()] = A_up;
        }
      }
      // free outflow boundary
      else {
        // TODO: make this more generic for other initial flow values
        const double A_init = param.A0;

        double W1, W2;

        if (in) {
          W1 = calculate_W1_value(Q_init, A_init, param.G0, param.rho, param.A0);
          W2 = calculate_W2_value(Q, A, param.G0, param.rho, param.A0);
        } else {
          W1 = calculate_W1_value(Q, A, param.G0, param.rho, param.A0);
          W2 = calculate_W2_value(Q_init, A_init, param.G0, param.rho, param.A0);
        }

        double Q_up = 0, A_up = 0;
        solve_W12(Q_up, A_up, W1, W2, param.G0, param.rho, param.A0);

        if (in) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_r[edge->get_id()] = A_up;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_l[edge->get_id()] = A_up;
        }
      }
    }
  }
}

void RightHandSideEvaluator::calculate_inner_fluxes(const std::vector<double> &u_prev) {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q_init = 0;

  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    if (vertex->get_edge_neighbors().size() == 2) {
      // get vertices
      const auto e0 = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      const auto e1 = d_graph->get_edge(vertex->get_edge_neighbors()[1]);

      // check orientation
      const bool e0_in = e0->is_pointing_to(vertex->get_id());
      const bool e1_in = e1->is_pointing_to(vertex->get_id());

      // get data
      const auto &data_e0 = e0->get_physical_data();
      const auto &data_e1 = e1->get_physical_data();
      const VesselParameters p_e0{data_e0.G0, data_e0.A0, data_e0.rho};
      const VesselParameters p_e1{data_e1.G0, data_e1.A0, data_e1.rho};

      // evaluate on edge e0
      const double Q_e0 =
        e0_in ? d_Q_macro_edge_boundary_value[2 * e0->get_id() + 1] : d_Q_macro_edge_boundary_value[2 * e0->get_id()];
      const double A_e0 =
        e0_in ? d_A_macro_edge_boundary_value[2 * e0->get_id() + 1] : d_A_macro_edge_boundary_value[2 * e0->get_id()];

      // evaluate on edge e1
      const double Q_e1 =
        e1_in ? d_Q_macro_edge_boundary_value[2 * e1->get_id() + 1] : d_Q_macro_edge_boundary_value[2 * e1->get_id()];
      const double A_e1 =
        e1_in ? d_A_macro_edge_boundary_value[2 * e1->get_id() + 1] : d_A_macro_edge_boundary_value[2 * e1->get_id()];

      assert(e0->get_vertex_neighbors()[1] == vertex->get_id() && e1->get_vertex_neighbors()[0] == vertex->get_id());

      const double W2_l = calculate_W2_value(Q_e0, A_e0, p_e0.G0, p_e0.rho, p_e0.A0);
      const double W1_r = calculate_W1_value(Q_e1, A_e1, p_e1.G0, p_e1.rho, p_e1.A0);

      assert(p_e0.G0 == p_e1.G0);
      assert(p_e0.rho == p_e1.rho);
      assert(p_e0.A0 == p_e1.A0);

      double Q_up = 0, A_up = 0;

      solve_W12(Q_up, A_up, W1_r, W2_l, p_e0.G0, p_e0.rho, p_e0.A0);

      d_Q_macro_edge_flux_r[e0->get_id()] = Q_up;
      d_Q_macro_edge_flux_l[e1->get_id()] = Q_up;
      d_A_macro_edge_flux_r[e0->get_id()] = A_up;
      d_A_macro_edge_flux_l[e1->get_id()] = A_up;
    }
  }
}

void RightHandSideEvaluator::calculate_fluxes(const double t, const std::vector<double> &u_prev) {
  calculate_nfurcation_fluxes(u_prev);
  // calculate_inner_fluxes(u_prev);
  calculate_inout_fluxes(t, u_prev);
}

void RightHandSideEvaluator::calculate_fluxes_on_macro_edge(const Edge &edge,
                                                            const std::vector<double> &u_prev,
                                                            std::vector<double> &Q_up_macro_edge,
                                                            std::vector<double> &A_up_macro_edge) {
  const auto local_dof_map = d_dof_map->get_local_dof_map(edge);

  assert(Q_up_macro_edge.size() == local_dof_map.num_micro_vertices());
  assert(A_up_macro_edge.size() == local_dof_map.num_micro_vertices());

  const auto &param = edge.get_physical_data();
  const double h = param.length / local_dof_map.num_micro_edges();

  const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

  // TODO: Make this more efficient by not recalculating the gradients.
  // finite-element for left and right edge
  FETypeNetwork fe(create_trapezoidal_rule(), num_basis_functions - 1);

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

  fe.reinit(h);

  // TODO: the boundary values of each cell are evaluated twice
  for (std::size_t micro_vertex_id = 1; micro_vertex_id < local_dof_map.num_micro_vertices() - 1; micro_vertex_id += 1) {
    const std::size_t local_micro_edge_id_l = micro_vertex_id - 1;
    const std::size_t local_micro_edge_id_r = micro_vertex_id;

    local_dof_map.dof_indices(local_micro_edge_id_l, 0, Q_dof_indices_l);
    local_dof_map.dof_indices(local_micro_edge_id_r, 0, Q_dof_indices_r);
    local_dof_map.dof_indices(local_micro_edge_id_l, 1, A_dof_indices_l);
    local_dof_map.dof_indices(local_micro_edge_id_r, 1, A_dof_indices_r);

    extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
    extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
    extract_dof(Q_dof_indices_r, u_prev, Q_prev_loc_r);
    extract_dof(A_dof_indices_r, u_prev, A_prev_loc_r);

    fe.evaluate_dof_at_quadrature_points(Q_prev_loc_l, Q_prev_qp_l);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc_l, A_prev_qp_l);
    fe.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

    const double Q_l = Q_prev_qp_l[1];
    const double A_l = A_prev_qp_l[1];
    const double Q_r = Q_prev_qp_r[0];
    const double A_r = A_prev_qp_r[0];

    const double W2_l = calculate_W2_value(Q_l, A_l, param.G0, param.rho, param.A0);
    const double W1_r = calculate_W1_value(Q_r, A_r, param.G0, param.rho, param.A0);

    double Q_up = 0, A_up = 0;

    solve_W12(Q_up, A_up, W1_r, W2_l, param.G0, param.rho, param.A0);

    Q_up_macro_edge[micro_vertex_id] = Q_up;
    A_up_macro_edge[micro_vertex_id] = A_up;
  }

  // update left fluxes
  Q_up_macro_edge[0] = d_Q_macro_edge_flux_l[edge.get_id()];
  A_up_macro_edge[0] = d_A_macro_edge_flux_l[edge.get_id()];

  // update right fluxes
  Q_up_macro_edge[local_dof_map.num_micro_vertices() - 1] = d_Q_macro_edge_flux_r[edge.get_id()];
  A_up_macro_edge[local_dof_map.num_micro_vertices() - 1] = d_A_macro_edge_flux_r[edge.get_id()];
}

void RightHandSideEvaluator::calculate_rhs(const double t, const std::vector<double> &u_prev, std::vector<double> &rhs) {
  // first zero rhs
  for (std::size_t idx = 0; idx < rhs.size(); idx += 1)
    rhs[idx] = 0;

  std::vector<double> Q_up_macro_edge(0, 0);
  std::vector<double> A_up_macro_edge(0, 0);

  // data structures for cell contributions
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto edge = d_graph->get_edge(e_id);
    const auto local_dof_map = d_dof_map->get_local_dof_map(*edge);

    // calculate fluxes on macro edge
    Q_up_macro_edge.resize(local_dof_map.num_micro_vertices());
    A_up_macro_edge.resize(local_dof_map.num_micro_vertices());
    calculate_fluxes_on_macro_edge(*edge, u_prev, Q_up_macro_edge, A_up_macro_edge);

    const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

    std::vector<double> f_loc_Q(num_basis_functions);
    std::vector<double> f_loc_A(num_basis_functions);

    QuadratureFormula qf = create_gauss4();
    FETypeNetwork fe(qf, local_dof_map.num_basis_functions() - 1);

    const auto &phi = fe.get_phi();
    const auto &phi_b = fe.get_phi_boundary();
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

    const auto &param = edge->get_physical_data();

    const double h = param.length / local_dof_map.num_micro_edges();

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      fe.reinit(h);
      qpm.reinit(micro_edge_id * h, static_cast<double>(1 + micro_edge_id) * h);

      // evaluate Q and A inside cell
      local_dof_map.dof_indices(micro_edge_id, 0, Q_dof_indices);
      local_dof_map.dof_indices(micro_edge_id, 1, A_dof_indices);

      extract_dof(Q_dof_indices, u_prev, Q_prev_loc);
      extract_dof(A_dof_indices, u_prev, A_prev_loc);

      fe.evaluate_dof_at_quadrature_points(Q_prev_loc, Q_prev_qp);
      fe.evaluate_dof_at_quadrature_points(A_prev_loc, A_prev_qp);

      // evaluate Q and A on boundary
      const double Q_up_0 = Q_up_macro_edge[micro_edge_id];
      const double Q_up_1 = Q_up_macro_edge[micro_edge_id + 1];
      const double A_up_0 = A_up_macro_edge[micro_edge_id];
      const double A_up_1 = A_up_macro_edge[micro_edge_id + 1];

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
        f_loc_Q[i] -= F_Q_eval(Q_up_1, A_up_1) * phi_b[1][i];
        f_loc_Q[i] += F_Q_eval(Q_up_0, A_up_0) * phi_b[0][i];

        f_loc_A[i] -= Q_up_1 * phi_b[1][i];
        f_loc_A[i] += Q_up_0 * phi_b[0][i];
      }

      // copy into global vector
      for (std::size_t i = 0; i < Q_dof_indices.size(); i += 1) {
        rhs[Q_dof_indices[i]] += f_loc_Q[i];
      }
      for (std::size_t i = 0; i < A_dof_indices.size(); i += 1) {
        rhs[A_dof_indices[i]] += f_loc_A[i];
      }
    }
  }
}

void RightHandSideEvaluator::apply_inverse_mass(std::vector<double> &rhs) {
  for (std::size_t i = 0; i < d_dof_map->num_dof(); i += 1)
    rhs[i] = d_inverse_mass[i] * rhs[i];
}

} // namespace macrocirculation