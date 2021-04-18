////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TRANSPORT_H
#define TUMORMODELS_TRANSPORT_H

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "flow_upwind_evaluator.hpp"
#include "gmm.h"
#include "graph_storage.hpp"
#include "right_hand_side_evaluator.hpp"

namespace macrocirculation {

class Transport {
public:
  Transport(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map_flow, std::shared_ptr<DofMap> dof_map_transport)
      : d_comm(comm),
        d_graph(std::move(graph)),
        d_dof_map_flow(std::move(dof_map_flow)),
        d_dof_map_transport(std::move(dof_map_transport)),
        d_flow_upwind_evaluator(comm, d_graph, d_dof_map_flow),
        d_edge_boundary_communicator(Communicator::create_edge_boundary_value_communicator(comm, d_graph)) {}

  void evaluate_macro_edge_boundary_values(const std::vector<double> &u_prev, const std::vector<double> &gamma_prev) {
    std::vector<std::size_t> dof_indices(4, 0);
    std::vector<double> local_dofs(4, 0);

    for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
      const auto edge = d_graph->get_edge(e_id);
      const auto &param = edge->get_physical_data();

      const auto &local_dof_map = d_dof_map_transport->get_local_dof_map(*edge);
      const double h = param.length / double(local_dof_map.num_micro_edges());

      FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
      fe.reinit(h);

      dof_indices.resize(local_dof_map.num_basis_functions());
      local_dofs.resize(local_dof_map.num_basis_functions());

      local_dof_map.dof_indices(0, 0, dof_indices);
      extract_dof(dof_indices, u_prev, local_dofs);
      d_gamma_macro_edge_boundary_value[2 * edge->get_id() + 0] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

      local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 0, dof_indices);
      extract_dof(dof_indices, u_prev, local_dofs);
      d_gamma_macro_edge_boundary_value[2 * edge->get_id() + 1] = fe.evaluate_dof_at_boundary_points(local_dofs).right;
    }
  }

  std::vector<double> &get_solution() { return d_solution; }

  void solve(double t, double dt, const std::vector<double> &u_prev) {
    evaluate_macro_edge_boundary_values(u_prev, d_solution);
    d_edge_boundary_communicator.update_ghost_layer(d_gamma_macro_edge_boundary_value);

    calculate_fluxes_at_nfurcations(t, u_prev);
    assemble_rhs(u_prev, d_rhs);
    apply_inverse_mass();
  }

  void calculate_fluxes_on_macro_edge(const double t,
                                      const Edge &edge,
                                      const std::vector<double> &u_prev,
                                      std::vector<double> &gamma_fluxes_edge) {

    const auto local_dof_map_flow = d_dof_map_flow->get_local_dof_map(edge);
    const auto local_dof_map_transport = d_dof_map_transport->get_local_dof_map(edge);

    assert(gamma_fluxes_edge.size() == local_dof_map_transport.num_micro_vertices());
    assert(gamma_fluxes_edge.size() == local_dof_map_flow.num_micro_vertices());
    // TODO: relax this precondition:
    assert(local_dof_map_flow.num_basis_functions() == local_dof_map_transport.num_micro_vertices());

    std::vector<double> Q_up_macro_edge(local_dof_map_flow.num_micro_vertices(), 0);
    std::vector<double> A_up_macro_edge(local_dof_map_flow.num_micro_vertices(), 0);

    d_flow_upwind_evaluator.get_fluxes_on_macro_edge(t, edge, u_prev, Q_up_macro_edge, A_up_macro_edge);

    const auto &param = edge.get_physical_data();
    const double h = param.length / double(local_dof_map_transport.num_micro_edges());

    const std::size_t num_basis_functions = local_dof_map_transport.num_basis_functions();

    // TODO: Make this more efficient by not recalculating the gradients.
    // finite-element for left and right edge
    FETypeNetwork fe(create_trapezoidal_rule(), num_basis_functions - 1);

    // dof indices for left and right edge
    std::vector<std::size_t> dof_indices(num_basis_functions, 0);
    std::vector<double> extracted_dof(num_basis_functions, 0);
    std::vector<double> evaluated_at_qps(num_basis_functions, 0);

    fe.reinit(h);

    // TODO: the boundary values of each cell are evaluated twice
    for (std::size_t micro_vertex_id = 1; micro_vertex_id < local_dof_map_flow.num_micro_vertices() - 1; micro_vertex_id += 1) {
      const std::size_t local_micro_edge_id_l = micro_vertex_id - 1;
      const std::size_t local_micro_edge_id_r = micro_vertex_id;

      // get gamma:
      local_dof_map_transport.dof_indices(local_micro_edge_id_l, 0, dof_indices);
      extract_dof(dof_indices, u_prev, extracted_dof);
      const double gamma_l = fe.evaluate_dof_at_boundary_points(extracted_dof).right;
      local_dof_map_transport.dof_indices(local_micro_edge_id_r, 0, dof_indices);
      extract_dof(dof_indices, u_prev, extracted_dof);
      const double gamma_r = fe.evaluate_dof_at_boundary_points(extracted_dof).left;

      const double v = Q_up_macro_edge[micro_vertex_id] / A_up_macro_edge[micro_vertex_id];
      const double gamma = (v < 0) ? gamma_r : gamma_l;

      gamma_fluxes_edge[micro_vertex_id] = gamma * v;
    }

    // update left fluxes
    gamma_fluxes_edge[0] = d_gamma_flux_l[edge.get_id()];

    // update right fluxes
    gamma_fluxes_edge[local_dof_map_transport.num_micro_vertices() - 1] = d_gamma_flux_r[edge.get_id()];
  }

  void assemble_rhs(double t, const std::vector<double> &u_prev, std::vector<double> &rhs) {
    // first zero rhs
    for (std::size_t idx = 0; idx < rhs.size(); idx += 1)
      rhs[idx] = 0;

    std::vector<double> gamma_flux_macro_edge(0, 0);

    // data structures for cell contributions
    for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
      const auto edge = d_graph->get_edge(e_id);
      const auto local_dof_map_flow = d_dof_map_flow->get_local_dof_map(*edge);
      const auto local_dof_map_transport = d_dof_map_transport->get_local_dof_map(*edge);

      // calculate fluxes on macro edge
      gamma_flux_macro_edge.resize(local_dof_map_transport.num_micro_vertices());
      calculate_fluxes_on_macro_edge(t, *edge, u_prev, gamma_flux_macro_edge);

      const std::size_t num_basis_functions = local_dof_map_transport.num_basis_functions();

      std::vector<double> f_loc_Q(num_basis_functions);
      std::vector<double> f_loc_A(num_basis_functions);

      QuadratureFormula qf = create_gauss4();
      FETypeNetwork fe(qf, local_dof_map_transport.num_basis_functions() - 1);

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

private:
  void calculate_fluxes_at_nfurcations(double t, const std::vector<double> &u_prev) {
    std::vector<double> Q_up_values(0, 0);
    std::vector<double> A_up_values(0, 0);

    for (auto &v_id : d_graph->get_vertex_ids()) {
      auto &vertex = *d_graph->get_vertex(v_id);

      d_flow_upwind_evaluator.get_fluxes_on_nfurcation(t, vertex, Q_up_values, A_up_values);

      if (vertex.is_leaf()) {
        auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);

        const double Q = Q_up_values.front();
        const double A = A_up_values.front();

        const double v = Q / A;

        const bool is_inflow = (v > 1e-8 && !edge.is_pointing_to(vertex.get_id())) || (v < 1e-8 && edge.is_pointing_to(vertex.get_id()));

        const double delta = 0.05;

        auto inflow_function = [=](double t) {
          if (vertex.is_inflow())
            return -2 * std::pow(t / delta, 3) + 3 * std::pow(t / delta, 2);
          else
            return 0.;
        };

        if (is_inflow) {
          // inflow
          if (!edge.is_pointing_to(vertex.get_id())) {
            d_gamma_flux_l[edge.get_id()] = Q * inflow_function(t);
          } else {
            d_gamma_flux_r[edge.get_id()] = Q * inflow_function(t);
          }
        } else {
          // outflow:
          if (!edge.is_pointing_to(vertex.get_id())) {
            d_gamma_flux_l[edge.get_id()] = v * d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 0];
          } else {
            d_gamma_flux_r[edge.get_id()] = v * d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 1];
          }
        }
      } else if (vertex.is_bifurcation()) {
        // do nothing if there is no flow
        if (gmm::vect_norminf(Q_up_values) < 1e-8)
          continue;

        std::vector<double> gamma;
        std::vector<bool> is_in;

        double Q_out = 0;
        double N_in = 0;

        for (size_t k = 0; k < vertex.get_edge_neighbors().size(); k += 1) {
          auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[k]);
          double gamma_value = edge.is_pointing_to(vertex.get_id()) ? d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 1] : d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 0];
          double v = Q_up_values[k] / A_up_values[k];
          const bool is_inflow_value = (v > 1e-8 && !edge.is_pointing_to(vertex.get_id())) || (v < 1e-8 && edge.is_pointing_to(vertex.get_id()));

          gamma.push_back(gamma_value);
          is_in.push_back((is_inflow_value));

          if (!is_inflow_value)
            Q_out += std::abs(Q_up_values[k]);

          if (is_inflow_value)
            N_in += std::abs(Q_up_values[k]) / A_up_values[k] * gamma_value;
        }

        for (std::size_t i = 0; i < vertex.get_edge_neighbors().size(); i += 1) {
          auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[i]);
          double v = Q_up_values[i] / A_up_values[i];

          if (is_in[i]) {
            if (edge.is_pointing_to(v_id))
              d_gamma_flux_r[edge.get_id()] = v * d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 1];
            else
              d_gamma_flux_l[edge.get_id()] = v * d_gamma_macro_edge_boundary_value[2 * edge.get_id() + 0];
          } else {
            double flux = v * N_in / Q_out * std::abs(Q_up_values[i]) * A_up_values[i] / std::abs(Q_up_values[i]);
            if (edge.is_pointing_to(v_id))
              d_gamma_flux_r[edge.get_id()] = flux;
            else
              d_gamma_flux_l[edge.get_id()] = flux;
          }
        }
      } else {
        throw std::runtime_error("not implemented");
      }
    }
  }

  void apply_inverse_mass();

private:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dof_map_flow;
  std::shared_ptr<DofMap> d_dof_map_transport;

  FlowUpwindEvaluator d_flow_upwind_evaluator;

  Communicator d_edge_boundary_communicator;

  std::vector<double> d_gamma_macro_edge_boundary_value;

  std::vector<double> d_gamma_flux_l;
  std::vector<double> d_gamma_flux_r;

  std::vector<double> d_rhs;
  std::vector<double> d_solution;
};

} // namespace macrocirculation

#endif //TUMORMODELS_TRANSPORT_H
