////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_GHOSTING_H
#define NETFV_GHOSTING_H

#include "utilLibs.hpp"

namespace netfv {

/**
 * Ghosting functor for 1-d elements
 */
class GhostingFunctorNet : public GhostingFunctor {
public:
  GhostingFunctorNet(const MeshBase &mesh)
      : d_mesh(mesh) {}

  virtual void operator()(const MeshBase::const_element_iterator &range_begin,
                          const MeshBase::const_element_iterator &range_end,
                          processor_id_type p, map_type &coupled_elements);

private:
  const MeshBase &d_mesh;
};

} // namespace netfv

#endif // NETFV_GHOSTING_H
