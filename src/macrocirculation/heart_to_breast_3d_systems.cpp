////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_3d_systems.hpp"

// makes the definition of HeartToBreast3DSolver precise
// HeartToBreast3DSolver is forward declared in include file heart_to_breast_3d_systems.hpp
#include "heart_to_breast_3d_solver.hpp"

namespace macrocirculation {

// define assembly functions
void CapillaryPressure::assemble() {}
void CapillaryPressure::assemble_1d() {}

void TissuePressure::assemble() {}

} // namespace macrocirculation