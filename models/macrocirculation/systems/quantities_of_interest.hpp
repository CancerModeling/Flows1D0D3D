////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_QUANTITIES_OF_INTEREST_HPP
#define TUMORMODELS_QUANTITIES_OF_INTEREST_HPP

#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class VesselDataStorage;
class DofMapNetwork;
template< std::size_t >
class FETypeInnerBdryNetwork;

template < std::size_t DEGREE >
void calculate_total_pressure(const GraphStorage &graph,
                             const VesselDataStorage& vessel_data,
                             const DofMapNetwork &map,
                             const FETypeInnerBdryNetwork<DEGREE> &fe,
                             const std::vector<double> &dof_vector,
                             std::vector<double> &interpolated);

}

#endif //TUMORMODELS_QUANTITIES_OF_INTEREST_HPP
