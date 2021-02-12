////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "vessel_data_storage.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

VesselDataStorage::VesselDataStorage()
    : d_parameters(0)
{}

std::size_t VesselDataStorage::add_parameter(const VesselParameters &parameters) {
  d_parameters.push_back(parameters);
  return d_parameters.size() - 1;
}

const VesselParameters &VesselDataStorage::get_parameters(const Edge &e) const {
  return d_parameters.at(e.get_type_id());
}

} // namespace macrocirculation
