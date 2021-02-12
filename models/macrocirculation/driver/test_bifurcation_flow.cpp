////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <cmath>

#include "../systems/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

inline double calculate_p_from_QA(double Q, double A, const mc::VesselParameters & param) {
  return 0.5*param.rho*std::pow(Q/A, 2) + param.G0 * ( std::sqrt(A/param.A0) - 1 );
}

void test1()
{
  const double Q_a = 1;
  const double Q_b = 0.5/2;
  const double Q_c = 0.5;

  const double A_a = 6.97;
  const double A_b = 6.97/2;
  const double A_c = 6.97;

  mc::VesselParameters param_a { 592.4e2,  6.97, 1.028 };
  mc::VesselParameters param_b { 592.4e2,  6.97/2, 1.028 };
  mc::VesselParameters param_c { 592.4e2,  6.97, 1.028 };

  double Q_a_up, A_a_up;
  double Q_b_up, A_b_up;
  double Q_c_up, A_c_up;

  const auto num_iter = mc::solve_at_bifurcation(
    Q_a, A_a, param_a, true,
    Q_b, A_b, param_b, false,
    Q_c, A_c, param_c, false,
    Q_a_up, A_a_up,
    Q_b_up, A_b_up,
    Q_c_up, A_c_up);

  std::cout << num_iter << " iterations" << std::endl;

  std::cout << "Q_a_up = " << Q_a_up << ", A_a_up = " << A_a_up << std::endl;
  std::cout << "Q_b_up = " << Q_b_up << ", A_b_up = " << A_b_up << std::endl;
  std::cout << "Q_c_up = " << Q_c_up << ", A_c_up = " << A_c_up << std::endl;

  std::cout << "diff Q = " << Q_a_up - Q_b_up - Q_c_up <<std::endl;
  std::cout << "diff p1 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_b_up, A_b_up, param_b) << std::endl;
  std::cout << "diff p2 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_c_up, A_c_up, param_c) << std::endl;
}

void test2()
{
  const double Q_a = 1;
  const double Q_b = 0.2;
  const double Q_c = 0.3;

  const double A_a = 7.2;
  const double A_b = 3.2;
  const double A_c = 4.1;

  mc::VesselParameters param_a { 592.4e2,  6.97, 1.028 };
  mc::VesselParameters param_b { 592.4e2,  3.1, 1.028 };
  mc::VesselParameters param_c { 592.4e2,  4.0, 1.028 };

  double Q_a_up, A_a_up;
  double Q_b_up, A_b_up;
  double Q_c_up, A_c_up;

  std::cout << "permutation true false false :" << std::endl;
  const auto num_iter = mc::solve_at_bifurcation(
    Q_a, A_a, param_a, true,
    Q_b, A_b, param_b, false,
    Q_c, A_c, param_c, false,
    Q_a_up, A_a_up,
    Q_b_up, A_b_up,
    Q_c_up, A_c_up);

  std::cout << num_iter << " iterations" << std::endl;

  std::cout << "Q_a_up = " << Q_a_up << ", A_a_up = " << A_a_up << std::endl;
  std::cout << "Q_b_up = " << Q_b_up << ", A_b_up = " << A_b_up << std::endl;
  std::cout << "Q_c_up = " << Q_c_up << ", A_c_up = " << A_c_up << std::endl;

  std::cout << "diff Q = " << Q_a_up - Q_b_up - Q_c_up <<std::endl;
  std::cout << "diff p1 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_b_up, A_b_up, param_b) << std::endl;
  std::cout << "diff p2 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_c_up, A_c_up, param_c) << std::endl;

  std::cout << "permutation false true false :" << std::endl;
  Q_a_up = A_a_up = Q_b_up = A_b_up = Q_c_up = A_c_up = 0;

  mc::solve_at_bifurcation(
    Q_b, A_b, param_b, false,
    Q_a, A_a, param_a, true,
    Q_c, A_c, param_c, false,
    Q_b_up, A_b_up,
    Q_a_up, A_a_up,
    Q_c_up, A_c_up);

  std::cout << "Q_a_up = " << Q_a_up << ", A_a_up = " << A_a_up << std::endl;
  std::cout << "Q_b_up = " << Q_b_up << ", A_b_up = " << A_b_up << std::endl;
  std::cout << "Q_c_up = " << Q_c_up << ", A_c_up = " << A_c_up << std::endl;

  std::cout << "diff Q = " << Q_a_up - Q_b_up - Q_c_up <<std::endl;
  std::cout << "diff p1 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_b_up, A_b_up, param_b) << std::endl;
  std::cout << "diff p2 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_c_up, A_c_up, param_c) << std::endl;

  std::cout << "permutation false false true:" << std::endl;
  Q_a_up = A_a_up = Q_b_up = A_b_up = Q_c_up = A_c_up = 0;

  mc::solve_at_bifurcation(
    Q_b, A_b, param_b, false,
    Q_c, A_c, param_c, false,
    Q_a, A_a, param_a, true,
    Q_b_up, A_b_up,
    Q_c_up, A_c_up,
    Q_a_up, A_a_up);

  std::cout << "Q_a_up = " << Q_a_up << ", A_a_up = " << A_a_up << std::endl;
  std::cout << "Q_b_up = " << Q_b_up << ", A_b_up = " << A_b_up << std::endl;
  std::cout << "Q_c_up = " << Q_c_up << ", A_c_up = " << A_c_up << std::endl;

  std::cout << "diff Q = " << Q_a_up - Q_b_up - Q_c_up <<std::endl;
  std::cout << "diff p1 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_b_up, A_b_up, param_b) << std::endl;
  std::cout << "diff p2 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_c_up, A_c_up, param_c) << std::endl;

  std::cout << "permutation false true true :" << std::endl;
  Q_a_up = A_a_up = Q_b_up = A_b_up = Q_c_up = A_c_up = 0;
  mc::solve_at_bifurcation(
    -Q_a, A_a, param_a, false,
    -Q_b, A_b, param_b, true,
    -Q_c, A_c, param_c, true,
    Q_a_up, A_a_up,
    Q_b_up, A_b_up,
    Q_c_up, A_c_up);

  std::cout << "Q_a_up = " << Q_a_up << ", A_a_up = " << A_a_up << std::endl;
  std::cout << "Q_b_up = " << Q_b_up << ", A_b_up = " << A_b_up << std::endl;
  std::cout << "Q_c_up = " << Q_c_up << ", A_c_up = " << A_c_up << std::endl;

  std::cout << "diff Q = " << Q_a_up - Q_b_up - Q_c_up <<std::endl;
  std::cout << "diff p1 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_b_up, A_b_up, param_b) << std::endl;
  std::cout << "diff p2 = " << calculate_p_from_QA(Q_a_up, A_a_up, param_a) - calculate_p_from_QA(Q_c_up, A_c_up, param_c) << std::endl;
}

int main(int argc, char *argv[]) {
  test1();
  test2();
}
