////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_UTIL_CREATE_3_VESSEL_NETWORK_H
#define TUMORMODELS_UTIL_CREATE_3_VESSEL_NETWORK_H

#include <memory>

// forward declarations;
namespace macrocirculation {
class GraphStorage;
}

namespace test_macrocirculation {

namespace util {

std::shared_ptr<macrocirculation::GraphStorage> create_3_vessel_network();

} // namespace util

} // namespace test_macrocirculation


#endif //TUMORMODELS_UTIL_CREATE_3_VESSEL_NETWORK_H
