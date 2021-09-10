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

#include "gmm.h"
#include "petsc.h"
#include <Eigen/Dense>

namespace macrocirculation {

class PetscMat {
public:
  PetscMat(const std::string &name, const std::vector<std::shared_ptr<DofMap>> &dof_map)
      : PetscMat(name,
                 std::accumulate(dof_map.begin(), dof_map.end(), 0, [](auto v, auto dm) { return v + dm->num_owned_dofs(); }),
                 std::accumulate(dof_map.begin(), dof_map.end(), 0, [](auto v, auto dm) { return v + dm->num_dof(); })) {}

  PetscMat(const std::string &name, const DofMap &dof_map)
      : PetscMat(name, dof_map.num_owned_dofs(), dof_map.num_dof()) {}

  PetscMat(const std::string &name, std::size_t num_owned_dofs, std::size_t num_dofs) {
    CHKERRABORT(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &d_mat));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetType(d_mat, MATMPIAIJ));
    CHKERRABORT(PETSC_COMM_WORLD, MatSetSizes(d_mat, num_owned_dofs, num_owned_dofs, num_dofs, num_dofs));
    // we overestimate the number non-zero entries, otherwise matrix assembly is incredibly slow :
    CHKERRABORT(PETSC_COMM_WORLD, MatMPIAIJSetPreallocation(d_mat, 1000, nullptr, 1000, nullptr));
    CHKERRABORT(PETSC_COMM_WORLD, PetscObjectSetName((PetscObject) d_mat, name.c_str()));
  }

  PetscMat(const PetscMat &) = delete;
  PetscMat(PetscMat &&) = delete;
  PetscMat &operator=(const PetscMat &) = delete;
  PetscMat &operator=(PetscMat &&) = delete;

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

  void add(const std::vector<size_t> &row_dofs, const std::vector<size_t> &col_dofs, Eigen::MatrixXd &mat_loc) {
    assert(row_dofs.size() == mat_loc.rows());
    assert(col_dofs.size() == mat_loc.cols());

    // we have no guarantees that gmm's memory is continuous
    // thus we copy it over
    std::vector<double> mat_memory(row_dofs.size() * col_dofs.size(), 0);
    for (int r = 0; r < row_dofs.size(); r += 1)
      for (int c = 0; c < col_dofs.size(); c += 1)
        mat_memory[c + r * col_dofs.size()] = mat_loc(r, c);

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

  double norm1() const {
    double value;
    CHKERRABORT(PETSC_COMM_WORLD, MatNorm(d_mat, NORM_1, &value));
    return value;
  }

  void zero() {
    MatZeroEntries(d_mat);
  }

  void assemble() {
    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyBegin(d_mat, MAT_FINAL_ASSEMBLY));
    CHKERRABORT(PETSC_COMM_WORLD, MatAssemblyEnd(d_mat, MAT_FINAL_ASSEMBLY));
  }

  Mat &get_mat() { return d_mat; }

  size_t first_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(PETSC_COMM_WORLD, MatGetOwnershipRange(d_mat, &start, &end));
    return start;
  }

  size_t last_dof() const {
    PetscInt start = 0;
    PetscInt end = 0;
    CHKERRABORT(PETSC_COMM_WORLD, MatGetOwnershipRange(d_mat, &start, &end));
    return end;
  }

  double get(size_t i, size_t j) const {
    double value;
    PetscInt ii = i;
    PetscInt jj = j;
    CHKERRABORT(PETSC_COMM_WORLD, MatGetValues(d_mat, 1, &ii, 1, &jj, &value));
    return value;
  }

private:
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
