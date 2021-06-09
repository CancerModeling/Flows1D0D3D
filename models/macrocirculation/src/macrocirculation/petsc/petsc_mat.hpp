////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_MAT_HPP
#define TUMORMODELS_PETSC_MAT_HPP

namespace macrocirculation {

class PetscMat {
public:
  // TODO: Rule of 5.

  PetscMat(const std::string &name, const DofMap &dof_map) {
    CHKERRABORT(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &d_mat));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetType(d_mat, MATMPIAIJ));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetType(d_mat, MATMPIAIJ));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetSizes(d_mat, (PetscInt) dof_map.num_owned_dofs(), dof_map.num_owned_dofs(), dof_map.num_dof(), dof_map.num_dof()));
    // we overestimate the number non-zero entries, otherwise matrix assembly is incredibly slow :
    CHKERRABORT(PETSC_COMM_WORLD, MatMPIAIJSetPreallocation(d_mat, 500, nullptr, 500, nullptr));
    // TODO: more generic name
    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_mat, name.c_str()));
  }

  void assemble(std::vector<PetscInt> rows, std::vector<PetscInt> cols, std::vector<double> values) {
    CHKERRABORT(PETSC_COMM_WORLD,
                MatSetValues(d_mat,
                             static_cast<PetscInt>(rows.size()),
                             rows.data(),
                             static_cast<PetscInt>(cols.size()),
                             cols.data(),
                             values.data(),
                             ADD_VALUES));

    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyBegin(d_mat, MAT_FINAL_ASSEMBLY));
    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyEnd(d_mat, MAT_FINAL_ASSEMBLY));
  }

  ~PetscMat() {
    CHKERRABORT(PETSC_COMM_WORLD, MatDestroy(&d_mat));
  }

  Mat &get_mat() { return d_mat; }

private:
  Mat d_mat;
};

}

#endif //TUMORMODELS_PETSC_MAT_HPP
