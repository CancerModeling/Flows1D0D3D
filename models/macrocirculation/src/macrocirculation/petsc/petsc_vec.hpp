////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_VEC_HPP
#define TUMORMODELS_PETSC_VEC_HPP

namespace macrocirculation {

class PetscVec {
public:
  // TODO: Rule of 5.

  PetscVec(const std::string &name, const DofMap &dof_map) {
    CHKERRABORT(PETSC_COMM_WORLD, VecCreate(PETSC_COMM_WORLD, &d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetType(d_vec, VECSTANDARD));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetSizes(d_vec, (PetscInt) dof_map.num_owned_dofs(), PETSC_DECIDE));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetUp(d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_vec, name.c_str()));
  }

  ~PetscVec() {
    CHKERRABORT(PETSC_COMM_WORLD, VecDestroy(&d_vec));
  }

  void set(PetscInt idx, double value) {
    CHKERRABORT(PETSC_COMM_WORLD, VecSetValue(d_vec, idx, value, INSERT_VALUES));
  }

  double get(PetscInt idx) {
    PetscReal value;
    CHKERRABORT(PETSC_COMM_WORLD, VecGetValues(d_vec, 1, &idx, &value));
    return value;
  }

  void assemble() {
    CHKERRABORT(PETSC_COMM_WORLD, VecAssemblyBegin(d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, VecAssemblyEnd(d_vec));
  }

private:
  Vec d_vec;
};

}

#endif //TUMORMODELS_PETSC_VEC_HPP
