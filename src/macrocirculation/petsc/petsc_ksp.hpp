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

  static std::shared_ptr<PetscKsp> create_named_from_options(PetscMat &mat, const std::string & name) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCJACOBI));

    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(ksp->d_ksp));

    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) ksp->d_ksp, (name).c_str()));

    return ksp;
  }

  static std::shared_ptr<PetscKsp> create_with_pc_ilu(PetscMat &mat) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    /*
    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCHYPRE));
    CHKERRABORT(PETSC_COMM_WORLD, PCHYPRESetType(pc, "euclid"));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetTolerances(ksp->d_ksp, 1e-16, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
     */

    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(ksp->d_ksp));
    // CHKERRABORT(PETSC_COMM_WORLD, KSPSetTolerances(ksp->d_ksp, 1e-16, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

    /*
    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc, PCJACOBI));
     */

    /*
    PC pc;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetType(ksp->d_ksp,KSPPREONLY));
    CHKERRABORT(PETSC_COMM_WORLD, PCSetType(pc,PCLU));
    CHKERRABORT(PETSC_COMM_WORLD, PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST));
     */
    // CHKERRABORT(PETSC_COMM_WORLD, KSPSetUp(ksp->d_ksp));

    return ksp;
  }

  void set_tolerances(double rtol, double atol, double divtol, int max_int)
  {
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetTolerances(d_ksp, rtol, atol, divtol, max_int));
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

  PetscInt get_iterations() const {
    PetscInt its;
    CHKERRABORT(PETSC_COMM_WORLD, KSPGetIterationNumber(d_ksp, &its));
    return its;
  }

protected:
  PetscKsp(PetscMat &mat) {
    CHKERRABORT(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &d_ksp));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetOperators(d_ksp, mat.get_mat(), mat.get_mat()));
    // CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_ksp, "pqsolver_"));
    CHKERRABORT(PETSC_COMM_WORLD, KSPSetFromOptions(d_ksp));
  }

private:
  KSP d_ksp;
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_KSP_HPP
