////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "petsc_assembly_blocks.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"

namespace macrocirculation {

Eigen::MatrixXd create_mass(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
  const auto &phi = fe.get_phi();
  const auto &JxW = fe.get_JxW();

  // TODO: This matrix is diagonal -> directly assemble it
  Eigen::MatrixXd m_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < static_cast<int>(local_dof_map.num_basis_functions()); j += 1) {
    for (int i = 0; i < static_cast<int>(local_dof_map.num_basis_functions()); i += 1) {
      m_loc(j, i) = 0;
      for (int qp = 0; qp < static_cast<int>(phi[i].size()); qp += 1)
        m_loc(j, i) += phi[i][qp] * phi[j][qp] * JxW[qp];
    }
  }
  return m_loc;
}

Eigen::MatrixXd create_phi_grad_psi(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  Eigen::MatrixXd k_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < static_cast<int>(local_dof_map.num_basis_functions()); j += 1) {
    for (int i = 0; i < static_cast<int>(local_dof_map.num_basis_functions()); i += 1) {
      k_loc(j, i) = 0;
      for (int qp = 0; qp < static_cast<int>(phi[i].size()); qp += 1)
        k_loc(j, i) += phi[i][qp] * dphi[j][qp] * JxW[qp];
    }
  }
  return k_loc;
}

Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row, BoundaryPointType col) {
  return create_boundary(local_dof_map.num_basis_functions(), row, col);
}

Eigen::MatrixXd create_boundary(size_t num_basis_functions, BoundaryPointType row, BoundaryPointType col) {
  Eigen::MatrixXd u_loc(num_basis_functions, num_basis_functions);
  const auto left = [](size_t i) -> double { return std::pow(-1., i); };
  const auto right = [](size_t) -> double { return 1.; };
  const auto phi = (col == BoundaryPointType::Left) ? left : right;
  const auto psi = (row == BoundaryPointType::Left) ? left : right;
  for (int j = 0; j < static_cast<int>(num_basis_functions); j += 1) {
    for (int i = 0; i < static_cast<int>(num_basis_functions); i += 1) {
      u_loc(j, i) = psi(j) * phi(i);
    }
  }
  return u_loc;
}

Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType type) {
  Eigen::VectorXd u_loc(local_dof_map.num_basis_functions());
  const auto left = [](size_t i) -> double { return std::pow(-1., i); };
  const auto right = [](size_t /*i*/) -> double { return 1.; };
  const auto phi = (type == BoundaryPointType::Left) ? left : right;
  for (int j = 0; j < static_cast<int>(local_dof_map.num_basis_functions()); j += 1)
    u_loc(j) = phi(j);
  return u_loc;
}

} // namespace macrocirculation