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

} // namespace macrocirculation

#endif // UTILS_H
