////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_ASSEMBLY_BLOCKS_HPP
#define TUMORMODELS_PETSC_ASSEMBLY_BLOCKS_HPP

#include "Eigen/Dense"

namespace macrocirculation {

// forward declarations:
class FETypeNetwork;
class LocalEdgeDofMap;

Eigen::MatrixXd create_mass(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map);

Eigen::MatrixXd create_phi_grad_psi(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map);

enum class BoundaryPointType { Left,
                               Right };

Eigen::MatrixXd create_boundary(size_t num_basis_functions, BoundaryPointType row, BoundaryPointType col);

Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row, BoundaryPointType col);

Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType type);

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_ASSEMBLY_BLOCKS_HPP
