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

default_S::default_S(double phi)
    : d_phi(phi) {}

void default_S::operator()(double,
                           const Edge &e,
                           const std::vector<double> &,
                           const std::vector<double> &Q,
                           const std::vector<double> &A,
                           std::vector<double> &S_Q_out,
                           std::vector<double> &S_A_out) const {
  // all vectors have to have the same shape
  assert(Q.size() == A.size());
  assert(Q.size() == S_Q_out.size());
  assert(Q.size() == S_A_out.size());

  const double mu = e.get_physical_data().viscosity;
  const double gamma = e.get_physical_data().gamma;

  for (std::size_t qp = 0; qp < Q.size(); qp += 1) {
    S_Q_out[qp] = -2 * mu * M_PI * (gamma + 2) * Q[qp] / A[qp];
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

  for (const auto &v_id : graph.get_active_vertex_ids(mpi::rank(comm))) {
    auto &vertex = *graph.get_vertex(v_id);

    // TODO: This is stupid!
    if (vertex.is_windkessel_outflow() || vertex.is_vessel_tree_outflow()) {
      auto &vertex_dof_map = dof_map.get_local_dof_map(vertex);
      auto &indices = vertex_dof_map.dof_indices();
      for (auto i : indices) {
        inv_mass[i] = 1;
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
      d_flow_upwind_evaluator(comm, d_graph, d_dof_map),
      d_S_evaluator(default_S{
        0 // 0 cm^2/s, no wall permeability
      }),
      d_degree(degree),
      d_inverse_mass(d_dof_map->num_dof()) {
  assemble_inverse_mass(d_comm, *d_graph, *d_dof_map, d_inverse_mass);
}

void RightHandSideEvaluator::evaluate(const double t, const std::vector<double> &u_prev, std::vector<double> &rhs) {
  d_flow_upwind_evaluator.init(t, u_prev);
  calculate_rhs(t, u_prev, rhs);
  // std::cout << "rhs " << rhs << std::endl;
  apply_inverse_mass(rhs);
  if (std::isnan(rhs.front()) || std::isnan(rhs.back())) {
    throw std::runtime_error("contains nans");
  }
}

void RightHandSideEvaluator::set_rhs_S(VectorEvaluator S_evaluator) {
  d_S_evaluator = std::move(S_evaluator);
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
    d_flow_upwind_evaluator.get_fluxes_on_macro_edge(t, *edge, u_prev, Q_up_macro_edge, A_up_macro_edge);

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
      d_S_evaluator(t, *edge, points, Q_prev_qp, A_prev_qp, S_Q, S_A);

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

  // add windkessel contributions
  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    Q_up_macro_edge.resize(1);
    A_up_macro_edge.resize(1);

    if (vertex->is_leaf() && vertex->is_windkessel_outflow()) {
      const auto &edge = *d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      assert(edge.has_physical_data());
      const auto &param = edge.get_physical_data();

      const bool is_pointing_to = edge.is_pointing_to(vertex->get_id());

      const auto &vertex_dof_map = d_dof_map->get_local_dof_map(*vertex);
      const auto &vertex_dofs = vertex_dof_map.dof_indices();

      d_flow_upwind_evaluator.get_fluxes_on_nfurcation(t, *vertex, Q_up_macro_edge, A_up_macro_edge);
      auto Q_out = Q_up_macro_edge.front();

      const double sgn = is_pointing_to ? +1 : -1;

      const auto p_c = u_prev[vertex_dofs[0]];

      // TODO: Move this calculation to the vertex.
      const double R1 = param.rho * param.get_c0() / param.A0;
      const double R2 = vertex->get_peripheral_vessel_data().resistance - R1;

      // pressure in the veins:
      const double p_v = vertex->get_peripheral_vessel_data().p_out;

      rhs[vertex_dofs[0]] = 1. / vertex->get_peripheral_vessel_data().compliance * (sgn * Q_out - (p_c - p_v) / R2);
    }
    // TODO: merge code with windkessel
    else if (vertex->is_leaf() && vertex->is_vessel_tree_outflow()) {
      const auto &edge = *d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      assert(edge.has_physical_data());

      const bool is_pointing_to = edge.is_pointing_to(vertex->get_id());

      const auto &vertex_dof_map = d_dof_map->get_local_dof_map(*vertex);
      const auto &vertex_dofs = vertex_dof_map.dof_indices();

      d_flow_upwind_evaluator.get_fluxes_on_nfurcation(t, *vertex, Q_up_macro_edge, A_up_macro_edge);
      auto Q_out = Q_up_macro_edge.front();

      const double sgn = is_pointing_to ? +1 : -1;

      std::vector<double> p_c;
      for (size_t k = 0; k < vertex_dofs.size(); k += 1) {
        p_c.push_back(u_prev[vertex_dofs[k]]);
      }

      // std::cout << p_c << std::endl;

      const auto &vtd = vertex->get_vessel_tree_data();
      assert(vertex_dofs.size() == vtd.capacitances.size());
      assert(vertex_dofs.size() == vtd.resistances.size());
      const auto &C = vtd.capacitances;
      // const double R_in = param.rho * param.get_c0() / param.A0;
      const double p_out = vtd.p_out;
      // we copy the resistances vector:
      std::vector<double> R{vtd.resistances};

      const auto n = static_cast< double >(vtd.furcation_number);

      // std::cout << "vid " << vertex->get_id() << " Q_out = " << Q_out << std::endl;

      // first
      rhs[vertex_dofs[0]] = 1. / C[0] * (sgn * Q_out - (p_c[0] - p_c[1]) / R[0]);

      // std::cout << rhs[vertex_dofs[0]] << std::endl;


      // middle
      for (size_t k = 1; k < p_c.size() - 1; k += 1)
        rhs[vertex_dofs[k]] = 1. / C[k] * ((p_c[k - 1] - p_c[k]) / (n * R[k - 1]) - (p_c[k] - p_c[k + 1]) / R[k]);

      //last
      const size_t k_last = p_c.size() - 1;
      rhs[vertex_dofs[k_last]] = 1. / C[k_last] * ((p_c[k_last - 1] - p_c[k_last]) / (n * R[k_last - 1]) - (p_c[k_last] - p_out) / (R[k_last]));
    }
  }
}

void RightHandSideEvaluator::apply_inverse_mass(std::vector<double> &rhs) {
  for (std::size_t i = 0; i < d_dof_map->num_dof(); i += 1)
    rhs[i] = d_inverse_mass[i] * rhs[i];
}

} // namespace macrocirculation