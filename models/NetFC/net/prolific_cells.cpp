////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "netUtil.hpp"
#include "../model.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

void netfc::Network::solve3DProlificCellProblem(int timeStep, double time) {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  // 3D ProlificCellProblem on a cube
  std::cout << " " << std::endl;
  std::cout << "Solve 3D ProlificCellProblem on a cube \Omega = (0," << L_x << ")^3" << std::endl;

  // time step size
  double dt = d_model_p->d_dt;
  std::cout << "dt: " << dt << std::endl;

  // Solver
  gmm::iteration iter(5.0e-11);

  // Switch solutions
  phi_P = phi_P_old;

  // Vectors
  std::vector<double> Delta_phi(2 * N_tot_3D, 0.0);

  // Newton iteration
  std::cout << " " << std::endl;
  std::cout << "Start Newton iteration 3D ProlificCellProblem" << std::endl;

  int N_newton = input.d_N_newton;

  std::cout << " " << std::endl;
  std::cout << "N_newton: " << N_newton << std::endl;

  double omega = input.d_omega;

  std::cout << " " << std::endl;
  std::cout << "omega: " << omega << std::endl;

  for (int i = 0; i < N_newton; i++) {

    std::cout << " " << std::endl;
    std::cout << "iteration: " << i << std::endl;
  }

  std::cout << " " << std::endl;
  std::cout << "Newton iteration 3D ProlificCellProblem finished" << std::endl;

  // Extract phi_p and mu_P
  for (int i = 0; i < N_tot_3D; i++) {

    phi_P[i] = sol_Vec_Phi_P[i];
    // mu_P[ i ] = sol_Vec_Phi_P[ i+N_tot_3D ];
  }

  // Switch solutions back
  phi_P_old = phi_P;
}
