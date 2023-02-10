////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2022 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <iostream>
#include <memory>

#include "macrocirculation/simple_linearized_solver.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "petsc.h"

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  {
    mc::SimpleLinearizedSolver solver_with_gap ("data/1d-meshes/vessels-with-gap.json", "vessels-with-gap" );
    mc::SimpleLinearizedSolver solver_gap ("data/1d-meshes/vessel-gap.json", "vessel-gap" );

    for (size_t i = 0; i < 1000000; i += 1) {
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
      if (i % 1000 == 0) {
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
