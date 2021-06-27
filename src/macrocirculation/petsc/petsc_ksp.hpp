////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "petsc_mat.hpp"
#include "petsc_vec.hpp"

#ifndef TUMORMODELS_PETSC_KSP_HPP
#define TUMORMODELS_PETSC_KSP_HPP

namespace macrocirculation {

class PetscKsp {
public:
  static std::shared_ptr<PetscKsp> create_with_pc_jacobi(PetscMat &mat) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCJACOBI));

    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(ksp->d_ksp));

    return ksp;
  }

  static std::shared_ptr<PetscKsp> create_with_pc_ilu(PetscMat &mat) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCHYPRE));
    CHKERRABORT(PETSC_COMM_WORLD, PCHYPRESetType(pc, "pilut"));

    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(ksp->d_ksp));

    return ksp;
  }

  PetscKsp(const PetscKsp &) = delete;
  PetscKsp(PetscKsp &&) = delete;
  PetscKsp &operator=(const PetscKsp &) = delete;
  PetscKsp &operator=(PetscKsp &&) = delete;

  ~PetscKsp() {
    CHKERRABORT(PETSC_COMM_WORLD, KSPDestroy(&d_ksp));
  }

  void solve(PetscVec &rhs, PetscVec &solution) {
    CHKERRABORT(PETSC_COMM_WORLD, KSPSolve(d_ksp, rhs.get_vec(), solution.get_vec()));
  }

  KSP &get_ksp() { return d_ksp; }

protected:
  PetscKsp(PetscMat &mat) {
    CHKERRABORT(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &d_ksp));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetOperators(d_ksp, mat.get_mat(), mat.get_mat()));
  }

private:
  KSP d_ksp;
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_KSP_HPP
