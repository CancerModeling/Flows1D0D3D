////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP
#define TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP

#include <vector>

namespace macrocirculation {

class GraphStorage;
class DofMapNetwork;

template < std::size_t DEGREE >
class FETypeNetwork;

template < std::size_t DEGREE >
void interpolate_to_vertices(const GraphStorage &graph,
                             const DofMapNetwork &map,
                             const FETypeNetwork<DEGREE> &fe,
                             const std::size_t component,
                             const std::vector<double> &dof_vector,
                             std::vector<double> &interpolated);

} // namespace macrocirculation

#endif //TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP
