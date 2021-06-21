////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include <cmath>

#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

inline double calculate_p_from_QA(double Q, double A, const mc::VesselParameters &param) {
  return 0.5 * param.rho * std::pow(Q / A, 2) + param.G0 * (std::sqrt(A / param.A0) - 1);
}

TEST_CASE("TestSimpleBifurcation", "[SimpleBifurcation]") {
  const double Q_a = 1;
  const double Q_b = 0.5 / 2;
  const double Q_c = 0.5;

  const double A_a = 6.97;
  const double A_b = 6.97 / 2;
  const double A_c = 6.97;

  mc::VesselParameters param_a{592.4e2, 6.97, 1.028};
  mc::VesselParameters param_b{592.4e2, 6.97 / 2, 1.028};
  mc::VesselParameters param_c{592.4e2, 6.97, 1.028};

  std::vector<double> Q_up(3, 0);
  std::vector<double> A_up(3, 0);

  const auto num_iter = mc::solve_at_nfurcation(
    {Q_a, Q_b, Q_c},
    {A_a, A_b, A_c},
    {param_a, param_b, param_c},
    {true, false, false},
    Q_up, A_up);

  INFO("needed " << num_iter << " iterations");

  INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
  INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
  INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);

  REQUIRE(Q_up[0] - Q_up[1] - Q_up[2] == Approx(0.).margin(1e-13));
  REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
  REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
}

TEST_CASE("TestPermutatedBifurcation", "[PermutatedBifurcation]") {
  const double Q_a = 1;
  const double Q_b = 0.2;
  const double Q_c = 0.3;

  const double A_a = 7.2;
  const double A_b = 3.2;
  const double A_c = 4.1;

  mc::VesselParameters param_a{592.4e2, 6.97, 1.028};
  mc::VesselParameters param_b{592.4e2, 3.1, 1.028};
  mc::VesselParameters param_c{592.4e2, 4.0, 1.028};

  {
    INFO("permutation true false false :");

    std::vector<double> Q_up(3, 0);
    std::vector<double> A_up(3, 0);

    const auto num_iter =
      mc::solve_at_nfurcation(
        {Q_a, Q_b, Q_c},
        {A_a, A_b, A_c},
        {param_a, param_b, param_c},
        {true, false, false},
        Q_up, A_up);

    INFO("needed " << num_iter << " iterations");

    INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
    INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
    INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);

    REQUIRE(Q_up[0] - Q_up[1] - Q_up[2] == Approx(0.).margin(1e-13));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
  }

  {
    INFO("permutation false true false :");

    std::vector<double> Q_up(3, 0);
    std::vector<double> A_up(3, 0);

    mc::solve_at_nfurcation(
      {Q_a, Q_b, Q_c},
      {A_a, A_b, A_c},
      {param_a, param_b, param_c},
      {false, true, false},
      Q_up, A_up);

    INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
    INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
    INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);

    REQUIRE(-Q_up[0] + Q_up[1] - Q_up[2] == Approx(0.).margin(1e-12));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
  }

  {
    INFO("permutation false false true:");

    std::vector<double> Q_up(3, 0);
    std::vector<double> A_up(3, 0);

    mc::solve_at_nfurcation(
      {Q_a, Q_b, Q_c},
      {A_a, A_b, A_c},
      {param_a, param_b, param_c},
      {false, false, true},
      Q_up, A_up);

    INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
    INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
    INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);

    REQUIRE(-Q_up[0] - Q_up[1] + Q_up[2] == Approx(0.).margin(1e-12));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
  }

  {
    INFO("permutation false true true :");

    std::vector<double> Q_up(3, 0);
    std::vector<double> A_up(3, 0);

    mc::solve_at_nfurcation(
      {Q_a, Q_b, Q_c},
      {A_a, A_b, A_c},
      {param_a, param_b, param_c},
      {false, true, true},
      Q_up, A_up);

    INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
    INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
    INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);

    REQUIRE(+Q_up[0] - Q_up[1] - Q_up[2] == Approx(0.).margin(1e-12));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
    REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
  }
}

TEST_CASE("TestTrifurcation", "[Trifurcation]") {
  const double Q_a = 1;
  const double Q_b = 0.5 / 2;
  const double Q_c = 0.5;
  const double Q_d = 0.1;

  const double A_a = 6.97;
  const double A_b = 6.97 / 2;
  const double A_c = 6.97;
  const double A_d = 6.97 / 5;

  mc::VesselParameters param_a{592.4e2, 6.97, 1.028};
  mc::VesselParameters param_b{592.4e2, 6.97 / 2, 1.028};
  mc::VesselParameters param_c{592.4e2, 6.97, 1.028};
  mc::VesselParameters param_d{592.4e2, 6.97 / 3, 1.028};

  std::vector<double> Q_up(4, 0);
  std::vector<double> A_up(4, 0);

  const auto num_iter =
    mc::solve_at_nfurcation(
      {Q_a, Q_b, Q_c, Q_d},
      {A_a, A_b, A_c, A_d},
      {param_a, param_b, param_c, param_d},
      {true, false, false, false},
      Q_up, A_up);

  INFO("needed " << num_iter << " iterations");

  INFO("Q_a_up = " << Q_up[0] << ", A_a_up = " << A_up[0]);
  INFO("Q_b_up = " << Q_up[1] << ", A_b_up = " << A_up[1]);
  INFO("Q_c_up = " << Q_up[2] << ", A_c_up = " << A_up[2]);
  INFO("Q_d_up = " << Q_up[3] << ", A_d_up = " << A_up[3]);

  REQUIRE(Q_up[0] - Q_up[1] - Q_up[2] - Q_up[3] == Approx(0.).margin(1e-12));

  REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[1], A_up[1], param_b) == Approx(0.).margin(1e-10));
  REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[2], A_up[2], param_c) == Approx(0.).margin(1e-10));
  REQUIRE(calculate_p_from_QA(Q_up[0], A_up[0], param_a) - calculate_p_from_QA(Q_up[3], A_up[3], param_d) == Approx(0.).margin(1e-10));
}
