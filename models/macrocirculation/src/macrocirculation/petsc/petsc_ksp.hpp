////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "petsc_mat.hpp"
#include "petsc_vec.hpp"

#ifndef TUMORMODELS_PETSC_KSP_HPP
#define TUMORMODELS_PETSC_KSP_HPP

namespace macrocirculation {

class PetscKsp {
public:
  PetscKsp(PetscMat &mat) {
    CHKERRABORT(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &d_ksp));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetOperators(d_ksp, mat.get_mat(), mat.get_mat()));

    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCJACOBI));

    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(d_ksp));
  }

  PetscKsp(const PetscKsp &) = delete;
  PetscKsp(PetscKsp &&) = delete;
  PetscKsp& operator=(const PetscKsp &) = delete;
  PetscKsp& operator=(PetscKsp &&) = delete;

  ~PetscKsp() {
    CHKERRABORT(PETSC_COMM_WORLD, KSPDestroy(&d_ksp));
  }

  void solve(PetscVec &rhs, PetscVec &solution) {
    CHKERRABORT(PETSC_COMM_WORLD, KSPSolve(d_ksp, rhs.get_vec(), solution.get_vec()));
  }

  KSP &get_ksp() { return d_ksp; }

private:
  KSP d_ksp;
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_KSP_HPP
