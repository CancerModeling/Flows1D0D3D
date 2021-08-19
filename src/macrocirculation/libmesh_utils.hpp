////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMESH_UTILS_H
#define LIBMESH_UTILS_H

#include "libmesh_includes.hpp"

namespace macrocirculation {

inline double get_min_nodal_spacing(const lm::MeshBase &mesh) {
  double h = 1.e20;
  for (const auto &nodei : mesh.local_node_ptr_range()) {
    for (const auto &nodej : mesh.local_node_ptr_range()) {
      auto dx = lm::Point(*nodei) - lm::Point(*nodej);
      auto d = dx.norm();
      if (d > 1.e-10 and d < h)
        h = d;
    }
  }

  mesh.comm().min(h);
  return h;
}

inline double get_mesh_size_estimate_using_element_volume(const lm::MeshBase &mesh) {
  double h = 1.e20;
  for (const auto &elem : mesh.active_element_ptr_range()) {
    double h2 = std::pow(3. * elem->volume() / (4. * M_PI), 1./3.);
    if (h2 < h)
      h = h2;
  }

  mesh.comm().min(h);
  return h;
}

} // namespace macrocirculation

#endif // LIBMESH_UTILS_H
