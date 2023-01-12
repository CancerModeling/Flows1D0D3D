////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_PETSC_MAT_HPP
#define TUMORMODELS_PETSC_MAT_HPP

#include <numeric>
#include <vector>

#include "../gmm_legacy_facade.hpp"
#include "petsc.h"
#include <Eigen/Dense>
#include "petsc_vec.hpp"

namespace macrocirculation {

class PetscMat {
public:
  PetscMat(MPI_Comm comm, const std::string &name, const std::vector<std::shared_ptr<DofMap>> &dof_map)
      : PetscMat(comm, 
                 name,
                 std::accumulate(dof_map.begin(), dof_map.end(), 0, [](auto v, auto dm) { return v + dm->num_owned_dofs(); }),
                 std::accumulate(dof_map.begin(), dof_map.end(), 0, [](auto v, auto dm) { return v + dm->num_dof(); })) {}

  PetscMat(MPI_Comm comm, const std::string &name, const DofMap &dof_map)
      : PetscMat(comm, name, dof_map.num_owned_dofs(), dof_map.num_dof()) {}

  PetscMat(MPI_Comm comm, const std::string &name, std::size_t num_owned_dofs, std::size_t num_dofs) : d_comm(comm) {
    CHKERRABORT(d_comm, MatCreate(d_comm, &d_mat));
    CHKERRABORT(d_comm, MatSetType(d_mat, MATMPIAIJ));
    CHKERRABORT(d_comm, MatSetSizes(d_mat, num_owned_dofs, num_owned_dofs, num_dofs, num_dofs));
    // we overestimate the number non-zero entries, otherwise matrix assembly is incredibly slow :
    CHKERRABORT(d_comm, MatMPIAIJSetPreallocation(d_mat, 1000, nullptr, 1000, nullptr));
    CHKERRABORT(d_comm, PetscObjectSetName((PetscObject) d_mat, name.c_str()));
  }

  PetscMat(const PetscMat &) = delete;
  PetscMat(PetscMat &&) = delete;
  PetscMat &operator=(const PetscMat &) = delete;
  PetscMat &operator=(PetscMat &&) = delete;

  ~PetscMat() {
    CHKERRABORT(d_comm, MatDestroy(&d_mat));
  }

  void add(const std::vector<size_t> &row_dofs, const std::vector<size_t> &col_dofs, gmm::row_matrix<gmm::wsvector<double>> &mat_loc) {
    assert(row_dofs.size() == mat_loc.nrows());
    assert(col_dofs.size() == mat_loc.ncols());

    // we have no guarantees that gmm's memory is continuous
    // thus we copy it over
    std::vector<double> mat_memory(row_dofs.size() * col_dofs.size(), 0);
    for (int r = 0; r < static_cast< int >( row_dofs.size() ); r += 1)
      for (int c = 0; c < static_cast< int >( col_dofs.size() ); c += 1)
        mat_memory[c + r * col_dofs.size()] = mat_loc[r][c];

    std::vector<PetscInt> row_dofs_(row_dofs.size(), 0);
    for (int r = 0; r < static_cast< int >( row_dofs.size() ); r += 1)
      row_dofs_[r] = static_cast<PetscInt>(row_dofs[r]);

    std::vector<PetscInt> col_dofs_(col_dofs.size(), 0);
    for (int c = 0; c < static_cast< int >( col_dofs.size() ); c += 1)
      col_dofs_[c] = static_cast<PetscInt>(col_dofs[c]);

    // pass it to petsc
    CHKERRABORT(d_comm,
                MatSetValues(d_mat,
                             static_cast<PetscInt>(row_dofs.size()),
                             row_dofs_.data(),
                             static_cast<PetscInt>(col_dofs.size()),
                             col_dofs_.data(),
                             mat_memory.data(),
                             ADD_VALUES));
  }

  void add(const std::vector<size_t> &row_dofs, const std::vector<size_t> &col_dofs, Eigen::MatrixXd &mat_loc) {
    assert(row_dofs.size() == mat_loc.rows());
    assert(col_dofs.size() == mat_loc.cols());

    // we have no guarantees that gmm's memory is continuous
    // thus we copy it over
    std::vector<double> mat_memory(row_dofs.size() * col_dofs.size(), 0);
    for (int r = 0; r < static_cast< int >( row_dofs.size() ); r += 1)
      for (int c = 0; c < static_cast< int >( col_dofs.size() ); c += 1)
        mat_memory[c + r * col_dofs.size()] = mat_loc(r, c);

    std::vector<PetscInt> row_dofs_(row_dofs.size(), 0);
    for (int r = 0; r < static_cast< int >( row_dofs.size() ); r += 1)
      row_dofs_[r] = static_cast<PetscInt>(row_dofs[r]);

    std::vector<PetscInt> col_dofs_(col_dofs.size(), 0);
    for (int c = 0; c < static_cast< int >( col_dofs.size() ); c += 1)
      col_dofs_[c] = static_cast<PetscInt>(col_dofs[c]);

    // pass it to petsc
    CHKERRABORT(d_comm,
                MatSetValues(d_mat,
                             static_cast<PetscInt>(row_dofs.size()),
                             row_dofs_.data(),
                             static_cast<PetscInt>(col_dofs.size()),
                             col_dofs_.data(),
                             mat_memory.data(),
                             ADD_VALUES));
  }

  double norm1() const {
    double value;
    CHKERRABORT(d_comm, MatNorm(d_mat, NORM_1, &value));
    return value;
  }

  void zero() {
    MatZeroEntries(d_mat);
  }

  void assemble() {
    CHKERRABORT(d_comm, MatAssemblyBegin(d_mat, MAT_FINAL_ASSEMBLY));
    CHKERRABORT(d_comm, MatAssemblyEnd(d_mat, MAT_FINAL_ASSEMBLY));
  }

  Mat &get_mat() { return d_mat; }

  size_t first_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(d_comm, MatGetOwnershipRange(d_mat, &start, &end));
    return start;
  }

  size_t last_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(d_comm, MatGetOwnershipRange(d_mat, &start, &end));
    return end;
  }

  double get(size_t i, size_t j) const {
    double value;
    PetscInt ii = i;
    PetscInt jj = j;
    CHKERRABORT(d_comm, MatGetValues(d_mat, 1, &ii, 1, &jj, &value));
    return value;
  }

  void mul(const PetscVec& src, PetscVec& dst) const {
    CHKERRABORT(d_comm, MatMult(d_mat, src.get_vec(), dst.get_vec()));
  }

  void print() const {
    CHKERRABORT(d_comm, MatView(d_mat, PETSC_VIEWER_STDOUT_WORLD));
    }

private:
  MPI_Comm d_comm;
  Mat d_mat{};
};

template<typename Stream>
inline Stream &operator<<(Stream &out, const PetscMat &m) {
  auto ldof = m.last_dof();
  for (size_t k = m.first_dof(); k < ldof; k += 1) {
    for (size_t j = m.first_dof(); j < ldof; j += 1)
      out << m.get(k, j) << " ";
    out << "\n";
  }
  return out;
}


} // namespace macrocirculation

#endif //TUMORMODELS_PETSC_MAT_HPP
