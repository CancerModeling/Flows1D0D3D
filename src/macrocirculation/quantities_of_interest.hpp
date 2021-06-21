////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_QUANTITIES_OF_INTEREST_HPP
#define TUMORMODELS_QUANTITIES_OF_INTEREST_HPP

#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class Point;

void calculate_total_pressure(MPI_Comm comm,
                              const GraphStorage &graph,
                              const DofMap &map,
                              const std::vector<double> &dof_vector,
                              std::vector<Point> &points,
                              std::vector<double> &interpolated);

void calculate_static_pressure(MPI_Comm comm,
                               const GraphStorage &graph,
                               const DofMap &map,
                               const std::vector<double> &dof_vector,
                               std::vector<Point> &points,
                               std::vector<double> &interpolated);

} // namespace macrocirculation

#endif //TUMORMODELS_QUANTITIES_OF_INTEREST_HPP
