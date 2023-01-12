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
    CHKERRABORT(mat.comm(), KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(mat.comm(), PCSetType(pc, PCJACOBI));

    CHKERRABORT(mat.comm(), KSPSetFromOptions(ksp->d_ksp));

    return ksp;
  }

  static std::shared_ptr<PetscKsp> create_named_from_options(PetscMat &mat, const std::string & name) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    PC pc;
    CHKERRABORT(mat.comm(), KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(mat.comm(), PCSetType(pc, PCJACOBI));

    CHKERRABORT(mat.comm(), KSPSetFromOptions(ksp->d_ksp));

    CHKERRABORT(mat.comm(), PetscObjectSetName((PetscObject) ksp->d_ksp, (name).c_str()));

    return ksp;
  }

  static std::shared_ptr<PetscKsp> create_with_pc_ilu(PetscMat &mat) {
    std::shared_ptr<PetscKsp> ksp(new PetscKsp(mat));

    /*
    PC pc;
    CHKERRABORT(d_comm, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(d_comm, PCSetType(pc, PCHYPRE));
    CHKERRABORT(d_comm, PCHYPRESetType(pc, "euclid"));
    CHKERRABORT(d_comm, KSPSetTolerances(ksp->d_ksp, 1e-16, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
     */

    PC pc;
    CHKERRABORT(mat.comm(), KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(mat.comm(), KSPSetFromOptions(ksp->d_ksp));
    // CHKERRABORT(d_comm, KSPSetTolerances(ksp->d_ksp, 1e-16, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

    /*
    PC pc;
    CHKERRABORT(d_comm, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(d_comm, PCSetType(pc, PCJACOBI));
     */

    /*
    PC pc;
    CHKERRABORT(d_comm, KSPGetPC(ksp->d_ksp, &pc));
    CHKERRABORT(d_comm, KSPSetType(ksp->d_ksp,KSPPREONLY));
    CHKERRABORT(d_comm, PCSetType(pc,PCLU));
    CHKERRABORT(d_comm, PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST));
     */
    // CHKERRABORT(d_comm, KSPSetUp(ksp->d_ksp));

    return ksp;
  }

  void set_tolerances(double rtol, double atol, double divtol, int max_int)
  {
    CHKERRABORT(d_comm, KSPSetTolerances(d_ksp, rtol, atol, divtol, max_int));
  }

  PetscKsp(const PetscKsp &) = delete;
  PetscKsp(PetscKsp &&) = delete;
  PetscKsp &operator=(const PetscKsp &) = delete;
  PetscKsp &operator=(PetscKsp &&) = delete;

  ~PetscKsp() {
    CHKERRABORT(d_comm, KSPDestroy(&d_ksp));
  }

  void solve(PetscVec &rhs, PetscVec &solution) {
    CHKERRABORT(d_comm, KSPSolve(d_ksp, rhs.get_vec(), solution.get_vec()));
  }

  KSP &get_ksp() { return d_ksp; }

  PetscInt get_iterations() const {
    PetscInt its;
    CHKERRABORT(d_comm, KSPGetIterationNumber(d_ksp, &its));
    return its;
  }

protected:
  PetscKsp(PetscMat &mat) : d_comm(mat.comm()) {
    CHKERRABORT(d_comm, KSPCreate(PETSC_COMM_WORLD, &d_ksp));
    CHKERRABORT(d_comm, KSPSetOperators(d_ksp, mat.get_mat(), mat.get_mat()));
    // CHKERRABORT(d_comm, PetscObjectSetName((PetscObject) d_ksp, "pqsolver_"));
    CHKERRABORT(d_comm, KSPSetFromOptions(d_ksp));
  }

private:
  MPI_Comm d_comm;
  KSP d_ksp;
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_KSP_HPP
