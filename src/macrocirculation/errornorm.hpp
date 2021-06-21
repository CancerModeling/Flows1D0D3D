////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_ERRORNORM_H
#define TUMORMODELS_ERRORNORM_H

#include <functional>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class FETypeNetwork;
class Point;

using EvaluationFunction = std::function<void(const std::vector<double> &, std::vector<double> &)>;

double errornorm(const MPI_Comm comm,
                 const GraphStorage &graph,
                 const DofMap &map,
                 std::size_t component,
                 const std::vector<double> &u_h,
                 const EvaluationFunction &u);

} // namespace macrocirculation

#endif //TUMORMODELS_ERRORNORM_H
