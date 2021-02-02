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
  explicit DofMapNetwork(std::size_t num_components)
      : d_num_components(num_components) {}

  void dof_indices(const mc::Edge &edge, std::vector<std::size_t> &dofs) const {
    const std::size_t num_basis_functions = 3; // since 2nd order
    dofs.resize(num_basis_functions * d_num_components);
    auto id = edge.get_id();
    for (auto idx = 0; idx < dofs.size(); idx += 1)
      dofs[idx] = id + idx;
  }

private:
  std::size_t d_num_components;
};

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  // create the ascending aorta
  mc::GraphStorage graph;
  std::size_t ascending_aorta_id = 1;
  auto v0 = graph.create_vertex(lm::Point(0, 0, 0));
  auto v1 = graph.create_vertex(lm::Point(4.0, 0, 0));
  graph.connect(*v0, *v1, ascending_aorta_id);
  graph.refine(4);

  const double velocity = 1;
  const double tau = 0.1;

  // assemble finite element system
  const std::size_t num_components = 1;
  const std::size_t num_basis_functions = 3;
  const std::size_t num_dofs = graph.num_edges() * num_components * num_basis_functions;

  gmm::row_matrix<gmm::wsvector<double>> A(num_dofs, num_dofs);
  std::vector<double> f(num_dofs);
  std::vector<double> u_now(num_dofs);
  std::vector<double> u_prev(num_dofs);

  DofMapNetwork dof_map(num_components);

  // reset matrix and rhs
  for (int i = 0; i < A.nrows(); i++)
    A[i].clear();

  for (int i = 0; i < f.size(); i++)
    f[i] = 0;

  // next time step
  std::swap(u_now, u_prev);

  // assemble cell integrals
  {
    gmm::row_matrix<gmm::wsvector<double>> A_loc(num_components * num_basis_functions, num_components * num_basis_functions);
    std::vector<double> f_loc(num_components * num_basis_functions);

    mc::FETypeNetwork fe(mc::create_gauss3());
    const auto &phi = fe.get_phi();
    const auto &dphi = fe.get_dphi();
    const auto &JxW = fe.get_JxW();

    std::vector<std::size_t> dof_indices;

    for (const auto &e_id : graph.get_edge_ids()) {
      const auto edge = graph.get_edge(e_id);
      fe.reinit(*edge);
      dof_map.dof_indices(*edge, dof_indices);

      // cell integral
      for (std::size_t i = 0; i < num_basis_functions; i += 1) {
        // rhs integral
        f_loc[i] = 0;
        for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
          f_loc[i] += phi[i][qp] * u_prev[i] * JxW[qp];
        }

        // matrix integral
        for (std::size_t j = 0; j < num_basis_functions; j += 1) {
          A_loc(i, j) = 0;
          for (std::size_t qp = 0; qp < phi[i].size(); qp += 1) {
            // mass for time derivative
            A_loc(i, j) += phi[i][qp] * phi[j][qp] * JxW[qp];
            // advection term
            A_loc(i, j) += (-tau) * velocity * phi[i][qp] * dphi[j][qp] * JxW[qp];
          }
        }
      }

      // copy into global matrix
      for (std::size_t i = 0; i < dof_indices.size(); i += 1) {
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

    // normals for inner boundaries
    const double normal_l = +1;
    const double normal_r = -1;

    // block matrices for inner boundaries
    gmm::row_matrix<gmm::wsvector<double>> A_ll_loc(num_components * num_basis_functions, num_components * num_basis_functions);
    gmm::row_matrix<gmm::wsvector<double>> A_lr_loc(num_components * num_basis_functions, num_components * num_basis_functions);
    gmm::row_matrix<gmm::wsvector<double>> A_rl_loc(num_components * num_basis_functions, num_components * num_basis_functions);
    gmm::row_matrix<gmm::wsvector<double>> A_rr_loc(num_components * num_basis_functions, num_components * num_basis_functions);

    // right hand side for inner boundaries
    std::vector<double> f_l_loc(num_components * num_basis_functions);
    std::vector<double> f_r_loc(num_components * num_basis_functions);

    // functions evaluated on the exterior boundaries
    mc::FETypeExteriorBdryNetwork fe_ext;
    const auto &phi = fe_ext.get_phi();

    // normal for the exterior boundaries
    const double normal_ext = -1;

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

        // zero local system
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          f_ext_loc[i] = 0;
          for (std::size_t j = 0; j < num_basis_functions; j += 1)
            A_ext_loc(i, j) = 0;
        }

        // inflow boundary
        if ((vertex->get_coordinate() - lm::Point(0, 0, 0)).norm() < 1e-14) {
        }
        // outflow boundary
        {
        }

        // copy into global matrix
        // TODO
      }
      // inner boundary
      else {
        const auto edge_l = graph.get_edge(vertex->get_edge_neighbors()[0]);
        const auto edge_r = graph.get_edge(vertex->get_edge_neighbors()[1]);
        dof_map.dof_indices(*edge_l, dof_indices_l);
        dof_map.dof_indices(*edge_r, dof_indices_r);

        // zero local system
        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          f_r_loc[i] = 0;
          f_l_loc[i] = 0;
          for (std::size_t j = 0; j < num_basis_functions; j += 1) {
            A_ll_loc(i, j) = 0;
            A_lr_loc(i, j) = 0;
            A_rl_loc(i, j) = 0;
            A_rr_loc(i, j) = 0;
          }
        }

        // ll
        for (std::size_t i = 0; i < num_basis_functions; i += 1)
          for (std::size_t j = 0; j < num_basis_functions; j += 1)
            if (false)
              true;

        // lr
        for (std::size_t i = 0; i < num_basis_functions; i += 1)
          for (std::size_t j = 0; j < num_basis_functions; j += 1)
            if (false)
              true;

        // rl
        for (std::size_t i = 0; i < num_basis_functions; i += 1)
          for (std::size_t j = 0; j < num_basis_functions; j += 1)
            if (false)
              true;

        // rr
        for (std::size_t i = 0; i < num_basis_functions; i += 1)
          for (std::size_t j = 0; j < num_basis_functions; j += 1)
            if (false)
              true;

        // copy into global matrix
        // TODO
      }
    }
  }

  return 0;
}
