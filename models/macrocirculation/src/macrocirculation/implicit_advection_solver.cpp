////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "implicit_advection_solver.hpp"

#include "graph_storage.hpp"

#include "dof_map.hpp"
#include "fe_type.hpp"
#include "gmm.h"
#include "graph_pvd_writer.hpp"
#include "interpolate_to_vertices.hpp"

namespace macrocirculation {

ImplicitAdvectionSolver::ImplicitAdvectionSolver(std::shared_ptr<GraphStorage> graph, std::size_t degree, double tau, double t_end, double velocity, InflowValueFct inflow_value_fct)
    : d_comm(MPI_COMM_WORLD),
      d_degree(degree),
      d_tau(tau),
      d_t_end(t_end),
      d_velocity(velocity),
      d_output_interval(1),
      d_inflow_value_fct(std::move(inflow_value_fct)),
      d_graph(std::move(graph)) {}

void ImplicitAdvectionSolver::solve() const {
  // assemble finite element system
  const std::size_t num_components = 1;
  const std::size_t num_basis_functions = d_degree + 1;
  const std::size_t num_dofs = d_graph->num_edges() * num_components * num_basis_functions;

  gmm::row_matrix<gmm::wsvector<double>> A(num_dofs, num_dofs);
  std::vector<double> f(num_dofs);
  std::vector<double> u_now(num_dofs, 0);
  std::vector<double> u_prev(num_dofs, 0);

  const auto dof_map = std::make_shared<DofMap>(d_graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *d_graph, num_components, d_degree, true);

  std::size_t it = 0;

  double t_now = 0;

  GraphPVDWriter writer(d_comm, "./output", "advection");

  while (t_now < d_t_end) {
    t_now += d_tau;
    it += 1;

    // reset matrix and rhs
    for (int i = 0; i < A.nrows(); i++)
      A[i].clear();

    for (int i = 0; i < f.size(); i++)
      f[i] = 0;

    // next time step
    for (int i = 0; i < u_now.size(); i++)
      u_prev[i] = u_now[i];

    {
      gmm::row_matrix<gmm::wsvector<double>> A_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      std::vector<double> f_loc(num_components * num_basis_functions);

      const auto qf = create_gauss4();

      FETypeNetwork fe(qf, d_degree);

      const auto &phi = fe.get_phi();
      const auto &dphi = fe.get_dphi();
      const auto &JxW = fe.get_JxW();

      std::vector<double> u_prev_loc(num_basis_functions, 0);
      std::vector<std::size_t> dof_indices(num_basis_functions);
      std::vector<double> u_prev_qp(qf.size(), 0);

      // assemble cell integrals
      for (const auto &e_id : d_graph->get_edge_ids()) {
        const auto macro_edge = d_graph->get_edge(e_id);

        const auto &local_dof_map = dof_map->get_local_dof_map(*macro_edge);
        const auto &param = macro_edge->get_physical_data();
        const double h = param.length / local_dof_map.num_micro_edges();

        fe.reinit(h);

        for (const auto &edge : macro_edge->micro_edges()) {
          local_dof_map.dof_indices(edge, 0, dof_indices);

          extract_dof(dof_indices, u_prev, u_prev_loc);
          fe.evaluate_dof_at_quadrature_points(u_prev_loc, u_prev_qp);

          // cell integral
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            // rhs integral
            f_loc[i] = 0;
            for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
              f_loc[i] += phi[i][qp] * u_prev_qp[qp] * JxW[qp];
            }

            // matrix integral
            for (std::size_t j = 0; j < num_basis_functions; j += 1) {
              A_loc(i, j) = 0;
              for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
                // mass for time derivative
                A_loc(i, j) += phi[i][qp] * phi[j][qp] * JxW[qp];
                // advection term
                A_loc(i, j) += (-d_tau) * d_velocity * phi[j][qp] * dphi[i][qp] * JxW[qp];
              }
            }
          }

          // copy into global matrix
          for (std::size_t i = 0; i < dof_indices.size(); i += 1) {
            f[dof_indices[i]] += f_loc[i];
            for (std::size_t j = 0; j < dof_indices.size(); j += 1) {
              A(dof_indices[i], dof_indices[j]) += A_loc(i, j);
            }
          }
        }
      }
    }

    // assemble inner boundary integrals
    {
      // functions evaluated on the inner cell boundaries
      FETypeNetwork fe(create_midpoint_rule(), d_degree);

      // phi_r is phi on the right edge seen from a given vertex
      const auto &phi_r = fe.get_phi_boundary()[0];
      // phi_l is phi on the left edge seen from a given vertex
      const auto &phi_l = fe.get_phi_boundary()[1];

      // block matrices for inner boundaries
      gmm::row_matrix<gmm::wsvector<double>> A_ll_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_lr_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_rl_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_rr_loc(num_components * num_basis_functions, num_components * num_basis_functions);

      // block matrices for exterior boundaries
      gmm::row_matrix<gmm::wsvector<double>> A_ext_loc(num_components * num_basis_functions, num_components * num_basis_functions);

      // right hand side for inner boundaries
      std::vector<double> f_ext_loc(num_components * num_basis_functions);

      std::vector<std::size_t> dof_indices_l;
      std::vector<std::size_t> dof_indices_r;

      for (const auto &e_id : d_graph->get_edge_ids()) {
        const auto macro_edge = d_graph->get_edge(e_id);

        const auto &local_dof_map = dof_map->get_local_dof_map(*macro_edge);
        const auto &param = macro_edge->get_physical_data();
        const double h = param.length / local_dof_map.num_micro_edges();

        fe.reinit(h);

        for (const auto &vertex : macro_edge->inner_micro_vertices()) {
          local_dof_map.dof_indices(*vertex.get_left_edge(), 0, dof_indices_l);
          local_dof_map.dof_indices(*vertex.get_right_edge(), 0, dof_indices_r);

          // zero local system
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            for (std::size_t j = 0; j < num_basis_functions; j += 1) {
              A_ll_loc(i, j) = 0;
              A_lr_loc(i, j) = 0;
              A_rl_loc(i, j) = 0;
              A_rr_loc(i, j) = 0;
            }
          }

          // ll
          if (d_velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_ll_loc(i, j) += d_tau * phi_l[j] * d_velocity * phi_l[i];
          }

          // lr
          if (d_velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_lr_loc(i, j) += d_tau * phi_r[j] * d_velocity * phi_l[i];
          }

          // rl
          if (d_velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rl_loc(i, j) += (-d_tau) * phi_l[j] * d_velocity * phi_r[i];
          }

          // rr
          if (d_velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rr_loc(i, j) += (-d_tau) * phi_r[j] * d_velocity * phi_r[i];
          }

          // copy into global matrix
          for (std::size_t i = 0; i < dof_indices_r.size(); i += 1) {
            for (std::size_t j = 0; j < dof_indices_r.size(); j += 1) {
              A(dof_indices_r[i], dof_indices_r[j]) += A_rr_loc(i, j);
              A(dof_indices_l[i], dof_indices_l[j]) += A_ll_loc(i, j);
              A(dof_indices_r[i], dof_indices_l[j]) += A_rl_loc(i, j);
              A(dof_indices_l[i], dof_indices_r[j]) += A_lr_loc(i, j);
            }
          }
        }

        for (const auto &edge : macro_edge->micro_edges()) {
        }
      }
    }

    // assemble outer boundaries
    {
    }


    /*
      for (const auto &v_id : d_graph->get_vertex_ids()) {
        const auto vertex = d_graph->get_vertex(v_id);

        // exterior boundary
        if (vertex->is_leaf()) {
          const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
          dof_map->dof_indices(*edge, dof_indices);
          fe_ext.reinit(*vertex, *edge);

          // zero local system
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            f_ext_loc[i] = 0;
            for (std::size_t j = 0; j < num_basis_functions; j += 1)
              A_ext_loc(i, j) = 0;
          }

          // inflow boundary
          if (fe_ext.get_normal() * d_velocity < 0) {
            const auto inflow_value = d_inflow_value_fct(t_now);
            for (unsigned int i = 0; i < num_basis_functions; i += 1) {
              f_ext_loc[i] += (-d_tau) * inflow_value * d_velocity * fe_ext.get_normal() * phi[i];
            }
          }
          // outflow boundary
          else {
            for (unsigned int i = 0; i < num_basis_functions; i += 1) {
              for (unsigned int j = 0; j < num_basis_functions; j += 1) {
                A_ext_loc(i, j) += d_tau * phi[i] * d_velocity * fe_ext.get_normal() * phi[j];
              }
            }
          }

          // copy into global matrix
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            f[dof_indices[i]] += f_ext_loc[i];
            for (std::size_t j = 0; j < num_basis_functions; j += 1) {
              A(dof_indices[i], dof_indices[j]) += A_ext_loc(i, j);
            }
          }
        }
        // inner boundary
        else {
          const auto edge_l = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
          const auto edge_r = d_graph->get_edge(vertex->get_edge_neighbors()[1]);
          dof_map->dof_indices(*edge_l, dof_indices_l);
          dof_map->dof_indices(*edge_r, dof_indices_r);
          fe_inner.reinit(*vertex, *edge_l, *edge_r);

          // zero local system
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            for (std::size_t j = 0; j < num_basis_functions; j += 1) {
              A_ll_loc(i, j) = 0;
              A_lr_loc(i, j) = 0;
              A_rl_loc(i, j) = 0;
              A_rr_loc(i, j) = 0;
            }
          }

          // ll
          if (d_velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_ll_loc(i, j) += d_tau * phi_l[j] * d_velocity * phi_l[i];
          }

          // lr
          if (d_velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_lr_loc(i, j) += d_tau * phi_r[j] * d_velocity * phi_l[i];
          }

          // rl
          if (d_velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rl_loc(i, j) += (-d_tau) * phi_l[j] * d_velocity * phi_r[i];
          }

          // rr
          if (d_velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rr_loc(i, j) += (-d_tau) * phi_r[j] * d_velocity * phi_r[i];
          }

          // copy into global matrix
          for (std::size_t i = 0; i < dof_indices.size(); i += 1) {
            for (std::size_t j = 0; j < dof_indices.size(); j += 1) {
              A(dof_indices_r[i], dof_indices_r[j]) += A_rr_loc(i, j);
              A(dof_indices_l[i], dof_indices_l[j]) += A_ll_loc(i, j);
              A(dof_indices_r[i], dof_indices_l[j]) += A_rl_loc(i, j);
              A(dof_indices_l[i], dof_indices_r[j]) += A_lr_loc(i, j);
            }
          }
        }
      }
    }
       */

    /*
    gmm::iteration iter(1.0E-10);
    gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A);
    gmm::gmres(A, u_now, f, PR, 500, iter);


    if (it % d_output_interval == 0) {
      std::cout << "it = " << it << std::endl;
      FETypeNetwork<DEGREE> fe(create_trapezoidal_rule());

      std::vector<double> u_vertex(d_graph->num_edges() * 2, 0);
      std::vector<Point> points(d_graph->num_edges() * 2);
      interpolate_to_vertices<DEGREE>(d_comm, *d_graph, *dof_map, 0, u_now, points, u_vertex);

      std::vector<double> u_mid(d_graph->num_edges(), 0);
      {
        for (std::size_t idx = 0; idx < d_graph->num_edges(); idx += 1) {
          u_mid[idx] = u_now[idx * num_basis_functions];
          if (num_basis_functions > 2)
            u_mid[idx] += -0.5 * u_now[idx * num_basis_functions + 2];
        }
      }

      writer.set_points(points);
      writer.add_vertex_data("concentration", u_vertex);
      writer.write(t_now);
    }
     */
  }
}

} // namespace macrocirculation