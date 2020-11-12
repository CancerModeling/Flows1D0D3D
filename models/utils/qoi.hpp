////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_QOI_H
#define UTILS_QOI_H

#include <iostream>
#include <vector>

namespace util {

struct QoI {

  std::string d_name;
  std::vector<double> d_data;

  QoI() {}
  QoI(const std::string &name) : d_name(name) {}

  void add(const double &data) { d_data.push_back(data); }

  const double &get_last() const {
    if (d_data.size() > 0)
      return d_data[d_data.size() - 1];
    else {

      libmesh_error_msg("Error: QoI data empty");
    }
  }

  double operator[](std::size_t i) const { return d_data[i]; }

  const std::vector<double> &get_all() const { return d_data; }
  std::vector<double> &get_all() { return d_data; }
};

struct QoIVec {

  std::vector<QoI> d_vec;

  QoIVec() {}

  QoIVec(const std::vector<std::string> &qoi_names) {
    for (const auto s : qoi_names)
      d_vec.emplace_back(s);
  }

  void add(const std::vector<double> &data) {
    for (unsigned int i = 0; i < data.size(); i++) {
      d_vec[i].add(data[i]);
    }
  }

  std::vector<std::string> get_names() const {

    std::vector<std::string> names;
    for (unsigned int i = 0; i < d_vec.size(); i++)
      names.emplace_back(d_vec[i].d_name);

    return names;
  }

  std::vector<double> get_last() const {
    auto data = std::vector<double>(d_vec.size(), 0.);
    for (unsigned int i = 0; i < d_vec.size(); i++) {
      data[i] = d_vec[i].get_last();
    }

    return data;
  }

  std::vector<std::vector<double>> get_all() const {
    std::vector<std::vector<double>> vec_data;

    for (unsigned int n = 0; n < d_vec[0].d_data.size(); n++) {

      auto data = std::vector<double>(d_vec.size(), 0.);

      for (unsigned int i = 0; i < d_vec.size(); i++) {
        data[i] = d_vec[i][n];
      }

      vec_data.emplace_back(data);
    }

    return vec_data;
  }
};

} // namespace util

#endif // UTILS_QOI_H
