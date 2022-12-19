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
#include "petsc.h"

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  {
    mc::SimpleLinearizedSolver solver;

    for (size_t i = 0; i < 10000; i += 1) {
      solver.solve();

      // output every 100
      if (i % 100 == 0) {
        // extract coupling data at aneurysm inflow
        auto in = solver.get_result(mc::SimpleLinearizedSolver::Outlet::in);
        std::cout << i << "  in: p = " << in.p << ", a = " << in.a << ", q = " << in.q << std::endl;

        // extract coupling data at aneurysm outflow
        auto out = solver.get_result(mc::SimpleLinearizedSolver::Outlet::out);
        std::cout << i << " out: p = " << out.p << ", a = " << out.a << ", q = " << out.q << std::endl;

        // just for fun, to see something, we could disable this
        solver.write();
      }
    }
  }

  CHKERRQ(PetscFinalize());
}
