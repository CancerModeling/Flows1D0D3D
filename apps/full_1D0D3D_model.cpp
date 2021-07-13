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

    // const double tau_out = tau;
    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    mc::HeartToBreast1DSolver solver(MPI_COMM_WORLD);
    solver.setup(degree, tau);

    // we average the pressure between the 9th and 10th heart beat:
    double t_start_pressure_averaging = static_cast<size_t>(std::floor(9. / tau));
    double t_stop_pressure_averaging = static_cast<size_t>(std::floor(10. / tau));

    const auto begin_t = std::chrono::steady_clock::now();
    double t = 0;
    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver.solve();

      if (it == t_start_pressure_averaging) {
        std::cout << "start integrating 0D pressures" << std::endl;
        solver.start_0d_pressure_integrator();
      }

      if (it == t_stop_pressure_averaging) {
        std::cout << "stop integrating 0D pressures" << std::endl;
        auto data = solver.stop_0d_pressure_integrator();

        for (auto &d : data) {

          // just return the values for now:
          std::cout << "rank = " << mc::mpi::rank(MPI_COMM_WORLD)
                    << ", x = " << d.p.x
                    << ", y = " << d.p.y
                    << ", z = " << d.p.z
                    << ", p_art = " << d.p_art
                    << ", p_ven = " << d.p_ven
                    << std::endl;
        }

        // Some condition to solve the 3D system
        {
          // TODO: Transfer 0D boundary values to 3D model.

          // TODO: Solver 3D system

          // TODO: Write 3D System
        }
      }

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", t = " << solver.get_time() << std::endl;

        solver.write_output();
      }

      // break
      if (t > t_end + 1e-12)
        break;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}
