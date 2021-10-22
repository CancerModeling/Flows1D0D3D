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

TEST_CASE("PressureAndFluxInflowBCAreConsistent", "[VesselFormulas]") {
  const double radius = 1.2;
  const double wall_thickness = 0.163;
  const double elastic_modulus = 400000.0;
  const double A0 = std::pow(radius, 2) * M_PI;

  const double G0 = mc::calculate_G0(wall_thickness, elastic_modulus, A0);
  const double rho = 1.028e-3;

  const double Q_DG = 10.;
  const double A_DG = 5.;

  const double Q_up = 12.;

  {
    const double A_up = mc::nonlinear::inflow::get_upwinded_A_from_Q(Q_DG, A_DG, true, Q_up, {G0, A0, rho});
    const double recreated_Q_up = mc::nonlinear::inflow::get_upwinded_Q_from_A(Q_DG, A_DG, +1., A_up, {G0, A0, rho});

    REQUIRE(recreated_Q_up == Approx(Q_up).epsilon(1e-8));
  }

  {
    const double A_up = mc::nonlinear::inflow::get_upwinded_A_from_Q(-Q_DG, A_DG, false, -Q_up, {G0, A0, rho});
    const double recreated_Q_up = mc::nonlinear::inflow::get_upwinded_Q_from_A(-Q_DG, A_DG, -1., A_up, {G0, A0, rho});

    REQUIRE(recreated_Q_up == Approx(-Q_up).epsilon(1e-8));
  }
}
