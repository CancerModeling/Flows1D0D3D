////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_BOUNDARY_CONDITION_READERS_HPP
#define TUMORMODELS_BOUNDARY_CONDITION_READERS_HPP

#include <string>

namespace macrocirculation {

// forward declarations:
class NonlinearLinearCoupling;

void read_coupling_conditions(NonlinearLinearCoupling & coupling, const std::string& filepath);

}

#endif //TUMORMODELS_BOUNDARY_CONDITION_READERS_HPP
