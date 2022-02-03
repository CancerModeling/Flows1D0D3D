////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include <utility>

#include "macrocirculation/upwinding_formulas_linearized_flow.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;


TEST_CASE("InnerBoundaryAndNfurcationBoundaryConincide", "[UpwindingFormulasLinearizedFlow]") {
  const double radius = 0.403;
  const double wall_thickness = 0.067;
  const double elastic_modulus = 400000.0;
  const double rho = 1.028e-3;

  const double A0 = std::pow(radius, 2) * M_PI;
  const double G0 = mc::calculate_G0(wall_thickness, elastic_modulus, A0);

  auto pdata = mc::VesselParameters(G0, A0, rho) ;

  const double C = mc::linear::get_C(pdata);
  const double L = mc::linear::get_L(pdata);

  double alpha = std::sqrt(C/L);

  double p_l = 100;
  double p_r = 105;
  double q_l = 0.15;
  double q_r = 0.10;
  double p_up = 0;
  double q_up = 0;

  mc::linearized::inner_boundary(alpha, p_l, q_l, p_r, q_r, p_up, q_up);

  const std::vector< double > p_vec = {p_l, p_r};
  const std::vector< double > q_vec = {q_l, q_r};
  std::vector< double > p_up_vec = {0, 0};
  std::vector< double > q_up_vec = {0, 0};

  // const std::vector< mc::VesselParameters> params = { {}, physical_data };

  mc::linearized::nfurcation_boundary(p_vec, q_vec, {pdata, pdata}, {+1, -1}, p_up_vec, q_up_vec);

  // test edges
  REQUIRE(p_up_vec[0] == Approx(p_up_vec[1]).epsilon(1e-12));
  REQUIRE(q_up_vec[0] == Approx(q_up_vec[1]).epsilon(1e-12));
  REQUIRE(p_up_vec[0] == Approx(p_up).epsilon(1e-12));
  REQUIRE(q_up_vec[0] == Approx(q_up).epsilon(1e-12));
}
