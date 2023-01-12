////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_VEC_HPP
#define TUMORMODELS_PETSC_VEC_HPP

#include <numeric>

#include "petsc.h"

namespace macrocirculation {

class PetscVec {
public:
  PetscVec(MPI_Comm comm, const std::string &name, const std::vector<std::shared_ptr<DofMap>> &dof_map)
      : PetscVec(comm, name, std::accumulate(dof_map.begin(), dof_map.end(), 0, [](auto v, auto dm) { return v + dm->num_owned_dofs(); })) {}

  PetscVec(MPI_Comm comm, const std::string &name, const DofMap &dof_map)
      : PetscVec(comm, name, dof_map.num_owned_dofs()) {}

  PetscVec(MPI_Comm comm, const std::string &name, size_t size)
      : d_comm(comm), d_vec() {
    CHKERRABORT(d_comm, VecCreate(d_comm, &d_vec));
    CHKERRABORT(d_comm, VecSetType(d_vec, VECSTANDARD));
    CHKERRABORT(d_comm, VecSetSizes(d_vec, (PetscInt) size, PETSC_DECIDE));
    CHKERRABORT(d_comm, VecSetUp(d_vec));
    CHKERRABORT(d_comm, PetscObjectSetName((PetscObject) d_vec, name.c_str()));
  }


  PetscVec(const PetscVec &) = delete;
  PetscVec(PetscVec &&) = delete;
  PetscVec &operator=(const PetscVec &) = delete;
  PetscVec &operator=(PetscVec &&) = delete;

  ~PetscVec() {
    CHKERRABORT(d_comm, VecDestroy(&d_vec));
  }

  void set_with_mode(const std::vector<size_t> &dofs, const std::vector<double> &values, InsertMode mode) {
    assert(dofs.size() == values.size());

    std::vector<PetscInt> dofs_(dofs.size(), 0);
    for (int r = 0; r < static_cast< int >( dofs.size() ); r += 1)
      dofs_[r] = static_cast<PetscInt>(dofs[r]);

    CHKERRABORT(d_comm, VecSetValues(d_vec, static_cast<PetscInt>(dofs_.size()), dofs_.data(), values.data(), mode));
  }

  void add(const std::vector<size_t> &dofs, const std::vector<double> &values) {
    set_with_mode(dofs, values, ADD_VALUES);
  }

  void set(const std::vector<size_t> &dofs, const std::vector<double> &values) {
    set_with_mode(dofs, values, INSERT_VALUES);
  }

  void set(PetscInt idx, double value) {
    CHKERRABORT(d_comm, VecSetValue(d_vec, idx, value, INSERT_VALUES));
  }

  double norm2() const {
    double value;
    CHKERRABORT(d_comm, VecNorm(d_vec, NORM_2, &value));
    return value;
  }

  double get(size_t idx) const { return get(static_cast<PetscInt>(idx)); }

  double get(PetscInt idx) const {
    PetscReal value;
    CHKERRABORT(d_comm, VecGetValues(d_vec, 1, &idx, &value));
    return value;
  }

  void assemble() {
    CHKERRABORT(d_comm, VecAssemblyBegin(d_vec));
    CHKERRABORT(d_comm, VecAssemblyEnd(d_vec));
  }

  void zero() {
    VecZeroEntries(d_vec);
  }

  Vec &get_vec() { return d_vec; }

  const Vec &get_vec() const { return d_vec; }

  void swap(PetscVec &v) {
    Vec tmp = v.d_vec;
    v.d_vec = d_vec;
    d_vec = tmp;
  }

  size_t local_size() const {
    PetscInt size = 0;
    CHKERRABORT(d_comm, VecGetLocalSize(d_vec, &size));
    return size;
  }

  size_t first_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(d_comm, VecGetOwnershipRange(d_vec, &start, &end));
    return start;
  }

  size_t last_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(d_comm, VecGetOwnershipRange(d_vec, &start, &end));
    return end;
  }

  void print() const {
    CHKERRABORT(d_comm, VecView(d_vec, PETSC_VIEWER_STDOUT_WORLD));
  }

private:
  MPI_Comm d_comm;
  Vec d_vec;
};

template<typename Stream>
inline Stream &operator<<(Stream &out, const PetscVec &v) {
  auto ldof = v.last_dof();
  for (size_t k = v.first_dof(); k < ldof; k += 1)
    out << v.get(k) << " ";
  return out;
}

} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_VEC_HPP
