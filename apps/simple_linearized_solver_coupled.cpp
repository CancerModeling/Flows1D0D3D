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
  return 1.1400127037858874 * sin(2 * M_PI * t) + 0.10000157404719523 * sin(4 * M_PI * t) - 0.11286043261149026 * sin(6 * M_PI * t) - 0.10272017307948093 * sin(8 * M_PI * t) - 0.0054169100548564046 * sin(10 * M_PI * t) - 0.02557859579598561 * sin(12 * M_PI * t) - 0.078047431190483588 * sin(14 * M_PI * t) - 0.022033974667631729 * sin(16 * M_PI * t) - 0.034738483688161952 * sin(18 * M_PI * t) - 0.045615688241977252 * sin(20 * M_PI * t) + 0.018999552643855375 * sin(22 * M_PI * t) - 0.0084982548911405227 * sin(24 * M_PI * t) + 0.0066704505306768112 * sin(26 * M_PI * t) + 0.0205896141914854 * sin(28 * M_PI * t) + 0.00056547376330811749 * sin(30 * M_PI * t) + 0.0079161507767753492 * sin(32 * M_PI * t) + 0.011957921955760809 * sin(34 * M_PI * t) + 0.0019131544051249399 * sin(36 * M_PI * t) + 0.0048081271514646045 * sin(38 * M_PI * t) - 0.48681323281506877 * cos(2 * M_PI * t) - 0.53309141008559258 * cos(4 * M_PI * t) - 0.20633570858479031 * cos(6 * M_PI * t) - 0.036670810452806915 * cos(8 * M_PI * t) - 0.014049767327029486 * cos(10 * M_PI * t) - 0.066663082737686313 * cos(12 * M_PI * t) - 0.015523473625300599 * cos(14 * M_PI * t) + 0.023855278352215264 * cos(16 * M_PI * t) - 0.016520131034039841 * cos(18 * M_PI * t) + 0.04615603374764362 * cos(20 * M_PI * t) + 0.031258225120331377 * cos(22 * M_PI * t) + 0.0052060140772440802 * cos(24 * M_PI * t) + 0.030280763300256783 * cos(26 * M_PI * t) + 0.00202669624626461 * cos(28 * M_PI * t) - 0.00026986726295702868 * cos(30 * M_PI * t) + 0.0091328960066964764 * cos(32 * M_PI * t) - 0.0025995007712682956 * cos(34 * M_PI * t) - 0.003740301110564225 * cos(36 * M_PI * t) + 0.0020915964477192118 * cos(38 * M_PI * t) + 3.0675015000000001;
}

class Interpolator {
public:
  void add(double ptime, double pvalue) {
    time.push_back(ptime);
    value.push_back(pvalue);
  }

  double operator()(double current_time) const {
    if (value.size() < 1)
      throw std::runtime_error("Not enough values to interpolate");

    double value_now = value[value.size() - 1];
    double value_prev = value[value.size() - 2];
    double time_now = time[time.size() - 1];
    double time_prev = time[time.size() - 2];
    double zeta = (current_time - time_prev) / (time_now - time_prev);
    double value_interp = value_prev + zeta * (value_now - value_prev);
    return value_interp;
  }

private:
  std::vector<double> time;
  std::vector<double> value;
};

class OutletInterpolator {
public:
  OutletInterpolator(int num_outlets)
      : value_interpolator(num_outlets) {}

  void add(int outlet, double time, double value) {
    assert(outlet < value_interpolator.size());
    value_interpolator[outlet].add(time, value);
  }

  double operator()(int outlet, double current_time) const {
    assert(outlet < value_interpolator.size());
    return value_interpolator[outlet](current_time);
  }

private:
  std::vector<Interpolator> value_interpolator;
};

int main(int argc, char *argv[]) {
  /// TODO: Sehr groesses E
  /// TODO: Geometrie beschreiben
  /// TODO: Fluesse immer in vessel stueck hinein

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  std::cout << "Communicator size: " << mc::mpi::size(MPI_COMM_WORLD) << std::endl;

  // Time step width during our swing in phase. Can be quiet coarse:
  const double tau_swing_in = 1e-1;
  // Time step width for the implicit solver with a gap:
  const double tau = 1e-3;
  // Time step width for the explicit solver solving the gap domain:
  const double tau_explicit = 2.5e-5;

  // Duration of our swing in phase :
  const double t_swing_in = 2;

  // Which domain type? 0 are two cylinders, 1 is a bifurcation, and 2 is a complex aneurysm geometry
  int domain_type = 2;

  // Choose the correct file paths for different domain types:
  std::string filepath_with_gap;
  std::string filepath_gap;
  switch (domain_type) {
    case 0:
      filepath_with_gap = "data/1d-meshes/vessels-with-gap.json";
      filepath_gap = "data/1d-meshes/vessel-gap.json";
      break;
    case 1:
      filepath_with_gap = "data/1d-meshes/bifurcation-with-gap.json";
      filepath_gap = "data/1d-meshes/bifurcation-gap.json";
      break;
    case 2:
      filepath_with_gap = "data/1d-meshes/Graph0-rcr-with-gap.json";
      filepath_gap = "data/1d-meshes/Graph0-rcr-gap.json";
      break;
    default:
      throw std::runtime_error("Unknown domain type.");
  }

  {
    // Solver for the outer network:
    mc::SimpleLinearizedSolver solver_with_gap(PETSC_COMM_WORLD, filepath_with_gap, "output", "vessels-with-gap", tau);
    // Solver for the inner network to achieve a swing in phase:
    mc::SimpleLinearizedSolver solver_gap_swing_in(PETSC_COMM_WORLD, filepath_gap, "output", "vessel-gap", tau);
    // Explicit solver for the inner network:
    mc::SimpleNonlinearSolver solver_gap(PETSC_COMM_WORLD, filepath_gap, "output", "vessel-gap", tau_explicit);

    solver_with_gap.set_inflow(inflow_aneurysm1);

    // Swing in phase (for this we just use the implicit solvers):
    for (size_t i = 0; i < int(t_swing_in / tau_swing_in); i += 1) {
      solver_with_gap.solve();
      solver_gap_swing_in.solve();
      const double t_now = (i + 1) * tau_swing_in;

      // couple the solver with a gap to the gap solver
      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        auto res = solver_with_gap.get_result(outlet_idx);
        solver_gap_swing_in.set_result(outlet_idx, res.p, res.q);
      }

      // couple the gap solver to the solver with a gap
      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        auto res = solver_gap_swing_in.get_result(outlet_idx);
        solver_with_gap.set_result(outlet_idx, res.p, res.q);
      }
    }

    // Reset the time
    solver_with_gap.set_t(0);

    // We have to couple an implicit ODE solver to an explict ODE solver. Hence, we need to interpolate between the values of the implicit solver.
    OutletInterpolator p_interpolator_with_gap(solver_gap.get_num_coupling_points());
    OutletInterpolator q_interpolator_with_gap(solver_gap.get_num_coupling_points());
    for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
      auto res = solver_with_gap.get_result(outlet_idx);
      p_interpolator_with_gap.add(outlet_idx, solver_with_gap.get_current_t(), res.p);
      q_interpolator_with_gap.add(outlet_idx, solver_with_gap.get_current_t(), res.q);
    }
    // We have to couple an explicit ODE solver to an implicit ODE solver. Hence, we have to extrapolate the values from the explicit solver.
    OutletInterpolator p_interpolator_gap(solver_gap.get_num_coupling_points());
    OutletInterpolator q_interpolator_gap(solver_gap.get_num_coupling_points());
    for (int outlet_idx = 0; outlet_idx < solver_gap.get_num_coupling_points(); outlet_idx += 1) {
      p_interpolator_gap.add(outlet_idx, 0., 0.);
      q_interpolator_gap.add(outlet_idx, 0., 0.);
    }

    // simulation phase:
    for (size_t i = 0; i < int(10 / tau); i += 1) {
      solver_with_gap.solve();

      // Put the results of the solver_with_gap into the interpolators:
      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        auto res = solver_with_gap.get_result(outlet_idx);
        p_interpolator_with_gap.add(outlet_idx, solver_with_gap.get_current_t(), res.p);
        q_interpolator_with_gap.add(outlet_idx, solver_with_gap.get_current_t(), res.q);
      }

      // Inner loop to solve with the explicit solver on a much smaller scale
      for (size_t j = 0; j < int(tau / tau_explicit); j += 1) {
        const double t_now = i * tau + j * tau_explicit;

        // Set interpolated coupling conditions, from the solver_with_gap to the solver_gap:
        for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1)
          solver_gap.set_result(outlet_idx, p_interpolator_with_gap(outlet_idx, t_now), q_interpolator_with_gap(outlet_idx, t_now));

        solver_gap.solve();

        // Put the results of the solver_gap into the interpolators:
        for (int outlet_idx = 0; outlet_idx < solver_gap.get_num_coupling_points(); outlet_idx += 1) {
          auto res = solver_gap.get_result(outlet_idx);
          p_interpolator_gap.add(outlet_idx, solver_gap.get_current_t(), res.p);
          q_interpolator_gap.add(outlet_idx, solver_gap.get_current_t(), res.q);
        }
      }

      // Couple to implicit solver: Set interpolated coupling conditions, from the solver_gap to the solver_with_gap:
      for (int outlet_idx = 0; outlet_idx < solver_with_gap.get_num_coupling_points(); outlet_idx += 1) {
        solver_with_gap.set_result(outlet_idx, p_interpolator_gap(outlet_idx, (i+2)*tau), q_interpolator_gap(outlet_idx, (i+2)*tau));
      }

      // Output not every time step 
      if ((i + 1) % int(1e-1 / tau) == 0) {
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
