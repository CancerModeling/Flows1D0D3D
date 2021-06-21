////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H
#define TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H

#include <memory>
#include <string>

namespace macrocirculation {

// forward declarations:
class GraphStorage;

void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph, const std::string &name_prefix);

} // namespace macrocirculation

#endif //TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H
