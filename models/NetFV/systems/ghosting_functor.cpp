////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "ghosting_functor.hpp"

void netfv::GhostingFunctorNet::operator()(
  const MeshBase::const_element_iterator &range_begin,
  const MeshBase::const_element_iterator &range_end, processor_id_type p,
  GhostingFunctor::map_type &coupled_elements) {

  CouplingMatrix *nullcm = nullptr;

  for (const auto &elem : as_range(range_begin, range_end)) {

    auto elem_id = elem->id();

    std::vector<unsigned int> elem_to_add;
    elem_to_add.push_back(elem_id);
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);
        elem_to_add.push_back(neighbor->id());
      }
    }

    for (auto e : elem_to_add) {
      auto n_elem = d_mesh.elem_ptr(e);
      if (n_elem->processor_id() != p)
        coupled_elements.emplace(std::make_pair(n_elem, nullcm));
    }
  }
}
