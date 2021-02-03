////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include <iostream>

#include "gmm.h"

#include "../systems/advection_assembly.hpp"
#include "../systems/fe_type_network.hpp"
#include "../systems/graph_storage.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

/// Simple hack to emulate a dof-map.
class DofMapNetwork {
public:
  explicit DofMapNetwork(std::size_t num_components, std::size_t num_basis_functions)
      : d_num_components(num_components), d_num_basis_functions(num_basis_functions) {}

  void dof_indices(const mc::Edge &edge, std::vector<std::size_t> &dofs) const {
    dofs.resize(d_num_basis_functions * d_num_components);
    auto id = edge.get_id();
    for (auto idx = 0; idx < dofs.size(); idx += 1)
      dofs[idx] = id * d_num_basis_functions * d_num_components + idx;
  }

private:
  std::size_t d_num_components;
  std::size_t d_num_basis_functions;
};

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  // create the ascending aorta
  mc::GraphStorage graph;

  std::size_t ascending_aorta_id = 1;

  auto v0 = graph.create_vertex(lm::Point(0, 0, 0));
  auto v1 = graph.create_vertex(lm::Point(1.0, 0, 0));
  graph.connect(*v0, *v1, ascending_aorta_id);
  graph.refine(5);

  // const std::size_t N = 32;
  // auto v_prev = graph.create_vertex(lm::Point(0, 0, 0));
  // for (std::size_t k = 0; k < N; k += 1) {
  //  auto v_next = graph.create_vertex(lm::Point((k + 1.) / N, 0, 0));
  //  graph.connect(*v_prev, *v_next, ascending_aorta_id);
  //  v_prev = v_next;
  // }

  const double velocity = 1;
  const double tau = 0.1;

  double t_now = 0;
  const double t_end = 1;

  auto inflow_boundary_value = [&t_now](auto &) { return std::sin(M_PI * t_now); };

  // assemble finite element system
  const std::size_t num_components = 1;
  const std::size_t num_basis_functions = 3;
  const std::size_t num_dofs = graph.num_edges() * num_components * num_basis_functions;

  gmm::row_matrix<gmm::wsvector<double>> A(num_dofs, num_dofs);
  std::vector<double> f(num_dofs);
  std::vector<double> u_now(num_dofs, 0);
  std::vector<double> u_prev(num_dofs, 0);

  DofMapNetwork dof_map(num_components, num_basis_functions);

  std::size_t it = 0;

  while (t_now < t_end) {
    t_now += tau;
    it += 1;

    // reset matrix and rhs
    for (int i = 0; i < A.nrows(); i++)
      A[i].clear();

    for (int i = 0; i < f.size(); i++)
      f[i] = 0;

    // next time step
    for (int i = 0; i < u_now.size(); i++)
      u_prev[i] = u_now[i];

    // assemble cell integrals
    {
      gmm::row_matrix<gmm::wsvector<double>> A_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      std::vector<double> f_loc(num_components * num_basis_functions);

      mc::FETypeNetwork fe(mc::create_gauss3());
      //mc::FETypeNetwork fe(mc::create_midpoint_rule());
      const auto &phi = fe.get_phi();
      const auto &dphi = fe.get_dphi();
      const auto &JxW = fe.get_JxW();

      std::vector<std::size_t> dof_indices;

      for (const auto &e_id : graph.get_edge_ids()) {
        const auto edge = graph.get_edge(e_id);
        fe.reinit(*edge);
        dof_map.dof_indices(*edge, dof_indices);

        std::vector<double> u_prev_loc(phi[0].size(), 0);
        for (std::size_t qp = 0; qp < phi[0].size(); qp += 1) {
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            u_prev_loc[qp] += phi[i][qp] * u_prev[dof_indices[i]];
          }
        }

        // cell integral
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          // rhs integral
          f_loc[i] = 0;
          for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
            f_loc[i] += phi[i][qp] * u_prev_loc[qp] * JxW[qp];
          }

          // matrix integral
          for (std::size_t j = 0; j < num_basis_functions; j += 1) {
            A_loc(i, j) = 0;
            for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
              // mass for time derivative
              A_loc(i, j) += phi[i][qp] * phi[j][qp] * JxW[qp];
              // advection term
              A_loc(i, j) += (-tau) * velocity * phi[j][qp] * dphi[i][qp] * JxW[qp];
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

    // assemble boundary integrals
    {
      // functions evaluated on the inner cell boundaries
      mc::FETypeInnerBdryNetwork fe_inner;
      const auto &phi_l = fe_inner.get_phi_l();
      const auto &phi_r = fe_inner.get_phi_r();

      // block matrices for inner boundaries
      gmm::row_matrix<gmm::wsvector<double>> A_ll_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_lr_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_rl_loc(num_components * num_basis_functions, num_components * num_basis_functions);
      gmm::row_matrix<gmm::wsvector<double>> A_rr_loc(num_components * num_basis_functions, num_components * num_basis_functions);

      // functions evaluated on the exterior boundaries
      mc::FETypeExteriorBdryNetwork fe_ext;
      const auto &phi = fe_ext.get_phi();

      // block matrices for exterior boundaries
      gmm::row_matrix<gmm::wsvector<double>> A_ext_loc(num_components * num_basis_functions, num_components * num_basis_functions);

      // right hand side for inner boundaries
      std::vector<double> f_ext_loc(num_components * num_basis_functions);

      std::vector<std::size_t> dof_indices;
      std::vector<std::size_t> dof_indices_l;
      std::vector<std::size_t> dof_indices_r;

      for (const auto &v_id : graph.get_vertex_ids()) {
        const auto vertex = graph.get_vertex(v_id);

        // exterior boundary
        if (vertex->is_leaf()) {
          const auto edge = graph.get_edge(vertex->get_edge_neighbors()[0]);
          dof_map.dof_indices(*edge, dof_indices);
          fe_ext.reinit(*vertex, *edge);

          // zero local system
          for (std::size_t i = 0; i < num_basis_functions; i += 1) {
            f_ext_loc[i] = 0;
            for (std::size_t j = 0; j < num_basis_functions; j += 1)
              A_ext_loc(i, j) = 0;
          }

          // inflow boundary
          if ((vertex->get_coordinate() - lm::Point(0, 0, 0)).norm() < 1e-14) {
            const auto inflow_value = inflow_boundary_value(vertex->get_coordinate());
            for (unsigned int i = 0; i < num_basis_functions; i += 1) {
              f_ext_loc[i] += (-tau) * inflow_value * velocity * fe_ext.get_normal() * phi[i];
            }
          }
          // outflow boundary
          else {
            for (unsigned int i = 0; i < num_basis_functions; i += 1) {
              for (unsigned int j = 0; j < num_basis_functions; j += 1) {
                A_ext_loc(i, j) += tau * phi[i] * velocity * fe_ext.get_normal() * phi[j];
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
          const auto edge_l = graph.get_edge(vertex->get_edge_neighbors()[0]);
          const auto edge_r = graph.get_edge(vertex->get_edge_neighbors()[1]);
          dof_map.dof_indices(*edge_l, dof_indices_l);
          dof_map.dof_indices(*edge_r, dof_indices_r);
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
          if (velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_ll_loc(i, j) += tau * phi_l[j] * velocity * phi_l[i];
          }

          // lr
          if (velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_lr_loc(i, j) += tau * phi_r[j] * velocity * phi_l[i];
          }

          // rl
          if (velocity * 1 > 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rl_loc(i, j) += (-tau) * phi_l[j] * velocity * phi_r[i];
          }

          // rr
          if (velocity * 1 < 0) {
            for (std::size_t i = 0; i < num_basis_functions; i += 1)
              for (std::size_t j = 0; j < num_basis_functions; j += 1)
                A_rr_loc(i, j) += (-tau) * phi_r[j] * velocity * phi_r[i];
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

    gmm::iteration iter(1.0E-10);
    gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A);
    gmm::gmres(A, u_now, f, PR, 500, iter);

    std::vector<double> u_mid(graph.num_edges(), 0);
    {
      for (std::size_t idx = 0; idx < graph.num_edges(); idx += 1) {
        u_mid[idx] = u_now[idx * num_basis_functions];
        if (num_basis_functions > 2)
          u_mid[idx] += -0.5 * u_now[idx * 3 + 2];
      }
    }

    mc::GraphDataWriter writer;
    writer.add_midpoint_data("concentration", u_mid);
    writer.write_vtk("concentration", graph, it);
  }

  return 0;
}
