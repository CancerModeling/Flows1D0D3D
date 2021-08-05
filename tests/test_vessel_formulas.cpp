////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include <utility>

#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

TEST_CASE("StaticPressureAndAreaAreInverse", "[VesselFormulas]") {
  double radius = 1.2;
  double wall_thickness = 0.163;
  double elastic_modulus = 400000.0;
  double A0 = std::pow(radius, 2) * M_PI;

  const double G0 = mc::calculate_G0(wall_thickness, elastic_modulus, A0);
  const double rho = 1.028e-3;

  const double A = 5.;

  const mc::VesselParameters params (G0, A0, rho);

  const double recreated_A = mc::nonlinear::get_A_from_p(params, mc::nonlinear::get_p_from_A(params, A));

  REQUIRE(recreated_A == Approx(A).epsilon(1e-16));
}
