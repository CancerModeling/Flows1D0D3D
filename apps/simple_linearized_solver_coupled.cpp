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

#include "macrocirculation/simple_linearized_solver.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "petsc.h"

namespace mc = macrocirculation;

double inflow(double t)
{
  // flow at the middle of vessel 22
  t = std::fmod(t, 1.);
  return 0.41543206934041998*sin(2*M_PI*t) + 0.011373623654789493*sin(4*M_PI*t) - 0.067330725324793395*sin(6*M_PI*t) - 0.04897078745454933*sin(8*M_PI*t) - 0.0018214247759830425*sin(10*M_PI*t) - 0.019937535008386593*sin(12*M_PI*t) - 0.01844597776677017*sin(14*M_PI*t) + 0.0011912928632729562*sin(16*M_PI*t) - 0.0082910209962541691*sin(18*M_PI*t) + 0.003781546492319121*sin(20*M_PI*t) + 0.0052424696925149764*sin(22*M_PI*t) + 0.0007945895297226625*sin(24*M_PI*t) - 0.2203282273590095*cos(2*M_PI*t) - 0.20258640381446483*cos(4*M_PI*t) - 0.085344073535552983*cos(6*M_PI*t) + 0.01217573129517773*cos(8*M_PI*t) - 0.001183996452239509*cos(10*M_PI*t) - 0.011310719833439547*cos(12*M_PI*t) + 0.013488225091287331*cos(14*M_PI*t) + 0.0071717162305028719*cos(16*M_PI*t) + 0.0056504458975141988*cos(18*M_PI*t) + 0.011203120977257584*cos(20*M_PI*t) + 0.0006885326606651587*cos(22*M_PI*t) + 0.0010044648362705904*cos(24*M_PI*t) + 1.1797162699999999;
}

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  std::cout << "size: " << mc::mpi::size(MPI_COMM_WORLD) << std::endl;

  {
    const double tau = 1e-3;

    mc::SimpleLinearizedSolver solver_with_gap ("data/1d-meshes/vessels-with-gap.json", "output", "vessels-with-gap", tau );
    mc::SimpleLinearizedSolver solver_gap ("data/1d-meshes/vessel-gap.json", "output", "vessel-gap", tau );

    solver_with_gap.set_inflow(inflow);
    //solver_with_gap.set_outflow_rcr(1.28e2, 5e-3);
    solver_with_gap.set_outflow_rcr(1e2, 1.75e-3);

    for (size_t i = 0; i < int(10/tau); i += 1) {
      // TODO: iterate
      solver_with_gap.solve();
      solver_gap.solve();

      /*
      {
        auto r_in = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
        auto r_out = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, r_in.p, r_in.q);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::in, r_out.p, r_out.q);
      }
       */

      {
        auto r_in = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
        solver_gap.set_result(mc::SimpleLinearizedSolver::Outlet::in, r_in.p, r_in.q);
        auto r_out = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
        solver_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, r_out.p, -r_out.q);
        // solver_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, 1., 1.);
      }

      {
        auto r_in = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::in, r_in.p, r_in.q);
        auto r_out = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, r_out.p, r_out.q);
        // solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, 1., 1.);
      }

      // output every 100
      if (i % int(1e-2/tau) == 0) {
        {
          // extract coupling data at aneurysm inflow
          auto in = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] " << i << " 1  in: p = " << in.p << ", a = " << in.a << ", q = " << in.q << std::endl;

          // extract coupling data at aneurysm outflow
          auto out = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] " << i << " 1 out: p = " << out.p << ", a = " << out.a << ", q = " << out.q << std::endl;
        }

        {
          // extract coupling data at aneurysm inflow
          auto in = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] " << i << " 2  in: p = " << in.p << ", a = " << in.a << ", q = " << in.q << std::endl;

          // extract coupling data at aneurysm outflow
          auto out = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
          std::cout << "[rank=" << mc::mpi::rank(MPI_COMM_WORLD) << "] " << i << " 2 out: p = " << out.p << ", a = " << out.a << ", q = " << out.q << std::endl;
        }

        // just for fun, to see something, we could disable this
        solver_with_gap.write();
        solver_gap.write();
      }
    }
  }

  CHKERRQ(PetscFinalize());
}
