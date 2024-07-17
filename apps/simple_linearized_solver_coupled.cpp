////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2022 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/simple_linearized_solver.hpp"
#include "macrocirculation/simple_nonlinear_solver.hpp"
#include "petsc.h"
#include <fstream>
#include <nlohmann/json.hpp>

namespace mc = macrocirculation;

double inflow(double t) {
  // flow at the middle of vessel 22
  t = std::fmod(t, 1.);
  return 0.41543206934041998 * sin(2 * M_PI * t) + 0.011373623654789493 * sin(4 * M_PI * t) - 0.067330725324793395 * sin(6 * M_PI * t) - 0.04897078745454933 * sin(8 * M_PI * t) - 0.0018214247759830425 * sin(10 * M_PI * t) - 0.019937535008386593 * sin(12 * M_PI * t) - 0.01844597776677017 * sin(14 * M_PI * t) + 0.0011912928632729562 * sin(16 * M_PI * t) - 0.0082910209962541691 * sin(18 * M_PI * t) + 0.003781546492319121 * sin(20 * M_PI * t) + 0.0052424696925149764 * sin(22 * M_PI * t) + 0.0007945895297226625 * sin(24 * M_PI * t) - 0.2203282273590095 * cos(2 * M_PI * t) - 0.20258640381446483 * cos(4 * M_PI * t) - 0.085344073535552983 * cos(6 * M_PI * t) + 0.01217573129517773 * cos(8 * M_PI * t) - 0.001183996452239509 * cos(10 * M_PI * t) - 0.011310719833439547 * cos(12 * M_PI * t) + 0.013488225091287331 * cos(14 * M_PI * t) + 0.0071717162305028719 * cos(16 * M_PI * t) + 0.0056504458975141988 * cos(18 * M_PI * t) + 0.011203120977257584 * cos(20 * M_PI * t) + 0.0006885326606651587 * cos(22 * M_PI * t) + 0.0010044648362705904 * cos(24 * M_PI * t) + 1.1797162699999999;
}

double inflow_aneurysm1(const double t) {
  return 1.1400127037858874*sin(2*M_PI*t) + 0.10000157404719523*sin(4*M_PI*t) - 0.11286043261149026*sin(6*M_PI*t) - 0.10272017307948093*sin(8*M_PI*t) - 0.0054169100548564046*sin(10*M_PI*t) - 0.02557859579598561*sin(12*M_PI*t) - 0.078047431190483588*sin(14*M_PI*t) - 0.022033974667631729*sin(16*M_PI*t) - 0.034738483688161952*sin(18*M_PI*t) - 0.045615688241977252*sin(20*M_PI*t) + 0.018999552643855375*sin(22*M_PI*t) - 0.0084982548911405227*sin(24*M_PI*t) + 0.0066704505306768112*sin(26*M_PI*t) + 0.0205896141914854*sin(28*M_PI*t) + 0.00056547376330811749*sin(30*M_PI*t) + 0.0079161507767753492*sin(32*M_PI*t) + 0.011957921955760809*sin(34*M_PI*t) + 0.0019131544051249399*sin(36*M_PI*t) + 0.0048081271514646045*sin(38*M_PI*t) - 0.48681323281506877*cos(2*M_PI*t) - 0.53309141008559258*cos(4*M_PI*t) - 0.20633570858479031*cos(6*M_PI*t) - 0.036670810452806915*cos(8*M_PI*t) - 0.014049767327029486*cos(10*M_PI*t) - 0.066663082737686313*cos(12*M_PI*t) - 0.015523473625300599*cos(14*M_PI*t) + 0.023855278352215264*cos(16*M_PI*t) - 0.016520131034039841*cos(18*M_PI*t) + 0.04615603374764362*cos(20*M_PI*t) + 0.031258225120331377*cos(22*M_PI*t) + 0.0052060140772440802*cos(24*M_PI*t) + 0.030280763300256783*cos(26*M_PI*t) + 0.00202669624626461*cos(28*M_PI*t) - 0.00026986726295702868*cos(30*M_PI*t) + 0.0091328960066964764*cos(32*M_PI*t) - 0.0025995007712682956*cos(34*M_PI*t) - 0.003740301110564225*cos(36*M_PI*t) + 0.0020915964477192118*cos(38*M_PI*t) + 3.0675015000000001;
}

int main(int argc, char *argv[]) {
  /// TODO: Sehr groesses E
  /// TODO: Geometrie beschreiben
  /// TODO: Fluesse immer in vessel stueck hinein

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  std::cout << "size: " << mc::mpi::size(MPI_COMM_WORLD) << std::endl;

  {
    // const double tau = 2.5e-6;
    const double tau = 2.5e-5;

    // mc::SimpleLinearizedSolver solver_with_gap ("data/1d-meshes/vessels-with-gap.json", "output", "vessels-with-gap", tau );
    // mc::SimpleLinearizedSolver solver_gap ("data/1d-meshes/vessel-gap.json", "output", "vessel-gap", tau );

    // mc::SimpleLinearizedSolver solver_with_gap(PETSC_COMM_SELF, "data/1d-meshes/bifurcation-with-gap.json", "output", "vessels-with-gap", tau);
    // mc::SimpleLinearizedSolver solver_gap(PETSC_COMM_SELF, "data/1d-meshes/bifurcation-gap.json", "output", "vessel-gap", tau);

    //mc::SimpleLinearizedSolver solver_with_gap(PETSC_COMM_SELF, "data/1d-meshes/Graph0-rcr-with-gap.json", "output", "vessels-with-gap", tau);
    //mc::SimpleLinearizedSolver solver_gap(PETSC_COMM_SELF, "data/1d-meshes/Graph0-rcr-gap.json", "output", "vessel-gap", tau);

    mc::SimpleLinearizedSolver solver_with_gap(PETSC_COMM_SELF, "data/1d-meshes/Graph0-rcr-with-gap.json", "output", "vessels-with-gap", tau);
    mc::SimpleNonlinearSolver solver_gap(PETSC_COMM_SELF, "data/1d-meshes/Graph0-rcr-gap.json", "output", "vessel-gap", tau);

    solver_with_gap.set_inflow(inflow_aneurysm1);
    //solver_with_gap.set_inflow(inflow_aneurysm1);

    for (size_t i = 0; i < int(10 / tau); i += 1) {
      // TODO: iterate
      solver_with_gap.solve();
      solver_gap.solve();
      const double t_now = (i + 1) * tau;

      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        auto res = solver_with_gap.get_result(outlet_idx);
        solver_gap.set_result(outlet_idx, res.p, res.q);
      }

      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        auto res = solver_gap.get_result(outlet_idx);
        solver_with_gap.set_result(outlet_idx, res.p, res.q);
      }

      // output every 100
      if ((i + 1) % int(1e-2 / tau) == 0) {
        for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
          // extract coupling data at aneurysm inflow
          auto in = solver_with_gap.get_result(outlet_idx);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] no-gap " << i << " " << outlet_idx << " in: p = " << in.p << ", a = " << in.a << ", q = " << in.q << std::endl;
        }

        for (int outlet_idx = 0; outlet_idx < solver_gap.get_num_coupling_points(); outlet_idx += 1) {
          // extract coupling data at aneurysm inflow
          auto in = solver_gap.get_result(outlet_idx);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] with gap " << i << " " << outlet_idx << "  in: p = " << in.p << ", a = " << in.a << ", q = " << in.q << std::endl;
        }

        // just for fun, to see something, we could disable this
        solver_with_gap.write();
        solver_gap.write();
      }
    }
  }

  CHKERRQ(PetscFinalize());
}
