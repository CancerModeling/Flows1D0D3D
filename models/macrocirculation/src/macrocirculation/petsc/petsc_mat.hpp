////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_MAT_HPP
#define TUMORMODELS_PETSC_MAT_HPP

#include <petsc.h>
#include <vector>

namespace macrocirculation {

class PetscMat {
public:
  PetscMat(const std::string &name, const DofMap &dof_map) {
    CHKERRABORT(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &d_mat));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetType(d_mat, MATMPIAIJ));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetSizes(d_mat, dof_map.num_owned_dofs(), dof_map.num_owned_dofs(), dof_map.num_dof(), dof_map.num_dof()));
    // we overestimate the number non-zero entries, otherwise matrix assembly is incredibly slow :
    CHKERRABORT(PETSC_COMM_WORLD, MatMPIAIJSetPreallocation(d_mat, 1000, nullptr, 1000, nullptr));
    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_mat, name.c_str()));
  }

  PetscMat(const PetscMat &) = delete;
  PetscMat(const PetscMat &&) = delete;
  void operator=(const PetscMat &) = delete;
  void operator=(const PetscMat &&) = delete;

  ~PetscMat() {
    CHKERRABORT(PETSC_COMM_WORLD, MatDestroy(&d_mat));
  }

  void add(const std::vector<size_t> &row_dofs, const std::vector<size_t> &col_dofs, gmm::row_matrix<gmm::wsvector<double>> &mat_loc) {
    assert(row_dofs.size() == mat_loc.nrows());
    assert(col_dofs.size() == mat_loc.ncols());

    // we have no guarantees that gmm's memory is continuous
    // thus we copy it over
    std::vector<double> mat_memory(row_dofs.size() * col_dofs.size(), 0);
    for (int r = 0; r < row_dofs.size(); r += 1)
      for (int c = 0; c < col_dofs.size(); c += 1)
        mat_memory[c + r * col_dofs.size()] = mat_loc[r][c];

    std::vector<PetscInt> row_dofs_(row_dofs.size(), 0);
    for (int r = 0; r < row_dofs.size(); r += 1)
      row_dofs_[r] = static_cast<PetscInt>(row_dofs[r]);

    std::vector<PetscInt> col_dofs_(col_dofs.size(), 0);
    for (int c = 0; c < col_dofs.size(); c += 1)
      col_dofs_[c] = static_cast<PetscInt>(col_dofs[c]);

    // pass it to petsc
    CHKERRABORT(PETSC_COMM_WORLD,
                MatSetValues(d_mat,
                             static_cast<PetscInt>(row_dofs.size()),
                             row_dofs_.data(),
                             static_cast<PetscInt>(col_dofs.size()),
                             col_dofs_.data(),
                             mat_memory.data(),
                             ADD_VALUES));
  }

  void zero() {
    MatZeroEntries(d_mat);
  }

  void assemble() {
    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyBegin(d_mat, MAT_FINAL_ASSEMBLY));
    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyEnd(d_mat, MAT_FINAL_ASSEMBLY));
  }

  Mat &get_mat() { return d_mat; }

private:
  Mat d_mat{};
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_MAT_HPP
