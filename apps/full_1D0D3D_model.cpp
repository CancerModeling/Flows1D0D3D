////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <macrocirculation/graph_csv_writer.hpp>
#include <utility>

#include "petsc.h"

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/coupled_explicit_implicit_1d_solver.hpp"
#include "macrocirculation/csv_vessel_tip_writer.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/heart_to_breast_1d_solver.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/nonlinear_linear_coupling.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                               //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.))) //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                     //
    ("tau-coup", "time step size for updating the coupling", cxxopts::value<double>()->default_value("1e-3"))         //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("10"))                      //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    const double t_end = args["t-end"].as<double>();
    const std::size_t max_iter = 160000000;

    const auto tau = args["tau"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();
    const auto tau_coup = args["tau-coup"].as<double>();

    // const double tau_out = tau;
    const auto output_interval = static_cast<std::size_t>(tau_out / tau);
    const auto coupling_interval = static_cast<std::size_t>(tau_coup / tau);

    mc::HeartToBreast1DSolver solver(MPI_COMM_WORLD);
    solver.setup(degree, tau, mc::BoundaryModel::DiscreteRCRTree);

    const auto begin_t = std::chrono::steady_clock::now();
    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver.solve();

      if (it % coupling_interval == 0) {
        std::cout << "calculates coupling " << std::endl;
        auto data = solver.get_vessel_tip_pressures();

        for (auto &d : data) {
          // just return the values for now:
          if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
            std::cout << d.p.x << ", " << d.p.y << ", " << d.p.z << ", " << d.pressure << ", " << d.R2 << std::endl;
        }

        // Some condition to solve the 3D system
        {
          // TODO: Transfer 0D boundary values to 3D model.

          // TODO: Solver 3D system

          // TODO: Write 3D System
        }

        // update the boundary conditions of the 1D system:
        {
          std::map<size_t, double> new_tip_pressures;
          for (auto &d : data) {
            // TODO: Replace this with something more meaningful.
            //       Note that 50 mmHg is much too much and just here to see the change from 30 mmHg which are the default.
            new_tip_pressures[d.vertex_id] = 50 * 1.3333;
          }
          solver.update_vessel_tip_pressures(new_tip_pressures);
        }
      }

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", t = " << solver.get_time() << std::endl;

        solver.write_output();
      }

      // break
      if (solver.get_time() > t_end + 1e-12)
        break;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}
