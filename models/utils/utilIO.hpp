////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_UTIL_IO_H
#define UTILS_UTIL_IO_H

#include <iostream>
#include <vector>

namespace util {

/*! @brief Provides geometrical methods such as point inside rectangle */
namespace io {

/*!
 * @brief Computes specified amount of tab
 */
inline std::string getTabS(int nt) {
  std::string tabS = "";
  for (int i = 0; i < nt; i++)
    tabS += "\t";

  return tabS;
}

inline std::string printStr(const bool flag, int nt = 0) {
  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS << (flag ? "true" : "false");
  return oss.str();
}

inline std::string printStr(const Point &x, int nt = 0) {
  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "(" << x(0) << ", " << x(1) << ", " << x(2) << ")";
  return oss.str();
}

template <class T>
inline std::string printStr(const std::vector<T> &list, int nt = 0) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << l;
    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
}

template <class T>
inline std::string printStr(const std::vector<std::vector<T>> &list,
                            int nt = 0) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << "(";
    for (size_t k = 0; k < l.size(); k++) {
      oss << l[k];
      if (k < l.size() - 1)
        oss << ", ";
    }
    oss << ")";

    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
}

template <>
inline std::string printStr(const std::vector<Point> &list, int nt) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << "(" << l(0) << ", " << l(1) << ", " << l(2) << ")";
    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
}



template <class T> inline void print(const std::vector<T> &list, int nt = 0) {

  std::cout << printStr(list, nt);
}

bool read_network_file(std::string filename, std::vector<Point> &nodes,
                       std::vector<std::vector<unsigned int>> &elems,
                       std::vector<unsigned int> &node_boundary_flag,
                       std::vector<bool> &node_branch_flag,
                       std::vector<double> &elem_radius,
                       std::vector<double> &elem_viscosity);
} // namespace io

} // namespace util

#endif // UTILS_UTIL_IO_H
