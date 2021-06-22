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
  PetscVec(const std::string &name, const DofMap &dof_map)
      : d_vec() {
    CHKERRABORT(PETSC_COMM_WORLD, VecCreate(PETSC_COMM_WORLD, &d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetType(d_vec, VECSTANDARD));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetSizes(d_vec, (PetscInt) dof_map.num_owned_dofs(), PETSC_DECIDE));
    CHKERRABORT(PETSC_COMM_WORLD, VecSetUp(d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_vec, name.c_str()));
  }

  PetscVec(const PetscVec &) = delete;
  PetscVec(PetscVec &&) = delete;
  PetscVec& operator=(const PetscVec &) = delete;
  PetscVec& operator=(PetscVec &&) = delete;

  ~PetscVec() {
    CHKERRABORT(PETSC_COMM_WORLD, VecDestroy(&d_vec));
  }

  void add(const std::vector<size_t> &dofs, const std::vector<double> &values) {
    assert(dofs.size() == values.size());

    std::vector<PetscInt> dofs_(dofs.size(), 0);
    for (int r = 0; r < dofs.size(); r += 1)
      dofs_[r] = static_cast<PetscInt>(dofs[r]);

    CHKERRABORT(PETSC_COMM_WORLD, VecSetValues(d_vec, static_cast<PetscInt>(dofs_.size()), dofs_.data(), values.data(), ADD_VALUES));
  }

  void set(PetscInt idx, double value) {
    CHKERRABORT(PETSC_COMM_WORLD, VecSetValue(d_vec, idx, value, INSERT_VALUES));
  }

  double get(size_t idx) const { return get(static_cast<PetscInt>(idx)); }

  double get(PetscInt idx) const {
    PetscReal value;
    CHKERRABORT(PETSC_COMM_WORLD, VecGetValues(d_vec, 1, &idx, &value));
    return value;
  }

  void assemble() {
    CHKERRABORT(PETSC_COMM_WORLD, VecAssemblyBegin(d_vec));
    CHKERRABORT(PETSC_COMM_WORLD, VecAssemblyEnd(d_vec));
  }

  void zero() {
    VecZeroEntries(d_vec);
  }

  Vec &get_vec() { return d_vec; }

  void swap(PetscVec &v) {
    Vec tmp = v.d_vec;
    v.d_vec = d_vec;
    d_vec = tmp;
  }

private:
  Vec d_vec;
};

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_VEC_HPP
