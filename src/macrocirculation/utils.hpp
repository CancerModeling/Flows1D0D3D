////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>
#include <algorithm>

namespace macrocirculation {

/*! @brief Compute time difference of two chrono objects */
inline float time_diff(std::chrono::steady_clock::time_point begin,
                std::chrono::steady_clock::time_point end) {

  return std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                               begin)
    .count();
}

/*! @brief Locate key in vector of keys */
template <class T>
inline long locate_in_set(const T &key, const std::vector<T> &set) {

  for (int i = 0; i < set.size(); i++)
    if (set[i] == key)
      return i;

  return -1;
}

/*! @brief Given 1d vector element id, compute 3d element id for image data */
inline std::vector<int> index_1d_3d(int I, std::vector<int> dim) {

  int i = I % dim[0];
  int a = std::floor(I / dim[0]);
  int j = a % dim[1];
  int k = std::floor(a / dim[1]);
  return {i, j, k};
}

/*! @brief Given 1d vector element id, compute 3d element id for image data */
inline void index_1d_3d(int I, std::vector<int> dim, std::vector<int> &I_3d) {
  I_3d.resize(3);
  I_3d[0] = I % dim[0];
  int a = std::floor(I / dim[0]);
  I_3d[1] = a % dim[1];
  I_3d[2] = std::floor(a / dim[1]);
}

/*! @brief Given 3d vector element id, compute 1d element id for image data */
inline int index_3d_1d(std::vector<int> I, std::vector<int> dim) {

  return I[0] + I[1] * dim[0] + I[2] * dim[0] * dim[1];
}

/*! @brief Given 3d vector element id, compute 1d element id for image data */
inline int index_3d_1d(int *I, int *dim) {

  return I[0] + I[1] * dim[0] + I[2] * dim[0] * dim[1];
}

/*! @brief Locates voxel (1d representation) for a given point */
inline int locate_voxel_1d(std::vector<double> x, std::vector<int> dim, double grid_size = 1.) {

  return int(x[0]/grid_size) + int(x[1]/grid_size) * dim[0] + int(x[2]/grid_size) * dim[0] * dim[1];
}

/*! @brief Locates voxel (3d representation) for a given point */
inline std::vector<int> locate_voxel_3d(std::vector<double> x, std::vector<int> dim, double grid_size = 1.) {

  return {int(x[0]/grid_size), int(x[1]/grid_size), int(x[2]/grid_size)};
}

/*! @brief print set in a list */
template <class T>
inline std::string print_str(const std::vector<T> &list) {
  std::ostringstream oss;
  oss << "(";
  for (size_t i=0; i<list.size(); i++) {
    oss << list[i];
    if (i < list.size() - 1)
      oss << ", ";
    else
      oss << ")\n";
  }
  return oss.str();
}

/*!
 * @brief sort based on list and return indices
 * see - https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
 */
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<class T>
inline T avg(const std::vector<T> &list) {

  if (list.size() == 0)
    return T(0);

  T avg = 0.;
  for (const auto &l : list)
    avg += l;
  avg = avg / (double(list.size()));

  return avg;
}

template<class T>
inline T std(const std::vector<T> &list) {

  if (list.size() == 0)
    return T(0);

  // get mean
  auto mean = avg(list);

  T dev = 0.;
  for (const auto &l : list)
    dev += (l - mean) * (l - mean);
  dev = dev / (T(list.size()));

  return std::sqrt(dev);
}

} // namespace macrocirculation

#endif // UTILS_H
