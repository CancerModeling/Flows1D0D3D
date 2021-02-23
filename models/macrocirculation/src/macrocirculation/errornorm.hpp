////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_ERRORNORM_H
#define TUMORMODELS_ERRORNORM_H

#include <vector>
#include <functional>
#include "libmesh/point.h"

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMapNetwork;
template<std::size_t degree>
class FETypeNetwork;

using libMesh::Point;

using EvaluationFunction = std::function<void(const std::vector<Point> &, std::vector<double> &)>;

template<std::size_t degree>
double errornorm(const GraphStorage &graph,
                 const DofMapNetwork &map,
                 std::size_t component,
                 const std::vector<double> &u_h,
                 EvaluationFunction u);

} // namespace macrocirculation

#endif //TUMORMODELS_ERRORNORM_H
