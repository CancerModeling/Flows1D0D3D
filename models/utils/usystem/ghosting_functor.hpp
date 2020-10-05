////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_GHOSTING_H
#define UTIL_GHOSTING_H

#include "utilLibs.hpp"

namespace util {

/**
 * Ghosting functor for 1-d elements
 */
class GhostingFunctorFV : public GhostingFunctor {
public:
  GhostingFunctorFV(const MeshBase &mesh)
      : d_mesh(mesh) {}

  virtual void operator()(const MeshBase::const_element_iterator &range_begin,
                          const MeshBase::const_element_iterator &range_end,
                          processor_id_type p, map_type &coupled_elements) {

    const CouplingMatrix *const null_mat = nullptr;
    for (const auto &elem : as_range(range_begin, range_end)) {
      coupled_elements.emplace(elem, null_mat);

      // add all elements to this element list (dense matrix)
      for (const auto &elem_remote : d_mesh.active_element_ptr_range())
        coupled_elements.emplace(elem_remote, null_mat);

      // add elements adjacent to this element
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr)
          coupled_elements.emplace(elem->neighbor_ptr(side), null_mat);
      }
    }
  }

private:
  const MeshBase &d_mesh;
};

} // namespace util

#endif // UTIL_GHOSTING_H
