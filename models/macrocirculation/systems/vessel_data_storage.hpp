////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_VESSEL_DATA_STORAGE_HPP
#define TUMORMODELS_VESSEL_DATA_STORAGE_HPP

#include "vessel_parameters.hpp"
#include <vector>

namespace macrocirculation {

// forward declarations
class Edge;

class VesselDataStorage {
public:
  VesselDataStorage();

  const VesselParameters &get_parameters(const Edge &e) const;

  std::size_t add_parameter(const VesselParameters &parameters);

private:
  std::vector<VesselParameters> d_parameters;
};

} // namespace macrocirculation

#endif //TUMORMODELS_VESSEL_DATA_STORAGE_HPP
