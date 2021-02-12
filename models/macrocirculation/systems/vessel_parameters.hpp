////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_VESSEL_PARAMETERS_HPP
#define TUMORMODELS_VESSEL_PARAMETERS_HPP

namespace macrocirculation {

/*! @brief Saves the material constants on a vessel. */
struct VesselParameters {
  VesselParameters() = default;

  VesselParameters(double G0, double A0, double rho)
      : G0(G0), A0(A0), rho(rho) {}

  double G0{};
  double A0{};
  double rho{};
};

} // namespace macrocirculation

#endif //TUMORMODELS_VESSEL_PARAMETERS_HPP
