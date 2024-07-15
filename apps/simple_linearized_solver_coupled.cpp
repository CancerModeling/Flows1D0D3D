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
#include <nlohmann/json.hpp>
#include <fstream>
#include "petsc.h"

namespace mc = macrocirculation;

double inflow(double t)
{
  // flow at the middle of vessel 22
  t = std::fmod(t, 1.);
  return 0.41543206934041998*sin(2*M_PI*t) + 0.011373623654789493*sin(4*M_PI*t) - 0.067330725324793395*sin(6*M_PI*t) - 0.04897078745454933*sin(8*M_PI*t) - 0.0018214247759830425*sin(10*M_PI*t) - 0.019937535008386593*sin(12*M_PI*t) - 0.01844597776677017*sin(14*M_PI*t) + 0.0011912928632729562*sin(16*M_PI*t) - 0.0082910209962541691*sin(18*M_PI*t) + 0.003781546492319121*sin(20*M_PI*t) + 0.0052424696925149764*sin(22*M_PI*t) + 0.0007945895297226625*sin(24*M_PI*t) - 0.2203282273590095*cos(2*M_PI*t) - 0.20258640381446483*cos(4*M_PI*t) - 0.085344073535552983*cos(6*M_PI*t) + 0.01217573129517773*cos(8*M_PI*t) - 0.001183996452239509*cos(10*M_PI*t) - 0.011310719833439547*cos(12*M_PI*t) + 0.013488225091287331*cos(14*M_PI*t) + 0.0071717162305028719*cos(16*M_PI*t) + 0.0056504458975141988*cos(18*M_PI*t) + 0.011203120977257584*cos(20*M_PI*t) + 0.0006885326606651587*cos(22*M_PI*t) + 0.0010044648362705904*cos(24*M_PI*t) + 1.1797162699999999;
}

int main(int argc, char *argv[]) {
  /// TODO: Sehr groesses E
  /// TODO: Geometrie beschreiben 
  /// TODO: Fluesse immer in vessel stueck hinein 

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  std::cout << std::ios::scientific;

  std::cout << "size: " << mc::mpi::size(MPI_COMM_WORLD) << std::endl;

  {
    const double tau = 1e-4;

    // mc::SimpleLinearizedSolver solver_with_gap ("data/1d-meshes/vessels-with-gap.json", "output", "vessels-with-gap", tau );
    // mc::SimpleLinearizedSolver solver_gap ("data/1d-meshes/vessel-gap.json", "output", "vessel-gap", tau );

    mc::SimpleLinearizedSolver solver_with_gap (PETSC_COMM_SELF, "data/1d-meshes/bifurcation-with-gap.json", "output", "vessels-with-gap", tau, true );
    mc::SimpleLinearizedSolver solver_gap (PETSC_COMM_SELF, "data/1d-meshes/bifurcation-gap.json", "output", "vessel-gap", tau, false );

    solver_with_gap.set_inflow(inflow);
    //solver_with_gap.set_outflow_rcr(1.28e2, 5e-3);
    // solver_with_gap.set_outflow_rcr(1e2, 1.75e-3);

    /*
    std::vector< double > t;
    std::vector< double > p_1_in;
    std::vector< double > p_1_out;
    std::vector< double > p_2_out;
    std::vector< double > p_2_in;
    std::vector< double > q_1_in;
    std::vector< double > q_1_out;
    std::vector< double > q_2_out;
    std::vector< double > q_2_in;
    std::vector< double > a_1_in;
    std::vector< double > a_1_out;
    std::vector< double > a_2_out;
    std::vector< double > a_2_in;
    std::vector< double > r_1_in;
    std::vector< double > r_1_out;
    std::vector< double > r_2_out;
    std::vector< double > r_2_in;
    std::vector< double > v_1_in;
    std::vector< double > v_1_out;
    std::vector< double > v_2_out;
    std::vector< double > v_2_in;
    */

    for (size_t i = 0; i < int(10/tau); i += 1) {
      // TODO: iterate
      solver_with_gap.solve();
      solver_gap.solve();
      const double t_now = (i+1)*tau;

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
        solver_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, r_out.p, r_out.q);
        // solver_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, 1., 1.);
      }

      {
        auto r_in = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::in, r_in.p, r_in.q);
        auto r_out = solver_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
        solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, r_out.p, r_out.q);
        // solver_with_gap.set_result(mc::SimpleLinearizedSolver::Outlet::out, 1., 1.);
      }

        // if (t_now >= 0. && (i+1) % int(1e-2/tau) == 0)
        {
          auto r1_in = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::in);
          auto r1_out = solver_with_gap.get_result_outer(mc::SimpleLinearizedSolver::Outlet::in);
          auto r2_out = solver_with_gap.get_result_outer(mc::SimpleLinearizedSolver::Outlet::out);
          auto r2_in = solver_with_gap.get_result(mc::SimpleLinearizedSolver::Outlet::out);
          
          /*

          t.push_back(t_now);
          p_1_in.push_back(r1_in.p); 
          p_1_out.push_back(r1_out.p); 
          p_2_out.push_back(r2_out.p );
          p_2_in.push_back(r2_in.p );
          q_1_in.push_back(r1_in.q); 
          q_1_out.push_back(r1_out.q); 
          q_2_out.push_back(-r2_out.q );
          q_2_in.push_back(-r2_in.q );
          a_1_in.push_back(r1_in.a); 
          a_1_out.push_back(r1_out.a); 
          a_2_out.push_back(r2_out.a );
          a_2_in.push_back(r2_in.a );
          r_1_in.push_back(std::sqrt(r1_in.a/M_PI)); 
          r_1_out.push_back(std::sqrt(r1_out.a/M_PI)); 
          r_2_out.push_back(std::sqrt(r2_out.a/M_PI ));
          r_2_in.push_back(std::sqrt(r2_in.a/M_PI ));
          v_1_in.push_back(r1_in.q/r1_in.a); 
          v_1_out.push_back(r1_out.q/r1_out.a); 
          v_2_out.push_back(-r2_out.q/r2_out.a );
          v_2_in.push_back(-r2_in.q/r2_in.a );

            std::ofstream f("coupling-values.json", std::ios::out);

            using json = nlohmann::json;

            json j = {
              {"version", "0.1"},
              { "vessel_data", {
                {"coupling_1_inner", {
                  {"id", "coupling_1_inner"},
                  {"gamma", 2},
                  {"x", solver_with_gap.get_points()[0]},
                  {"p", p_1_in},
                  {"q", q_1_in},
                  {"r", r_1_in},
                  {"v", v_1_in},
                  {"a", a_1_in}}},
                {"coupling_1_outer", {
                  {"id", "coupling_1_outer"},
                  {"gamma", 2},
                  {"x", solver_with_gap.get_points()[1]},
                  {"p", p_1_out},
                  {"q", q_1_out},
                  {"r", r_1_out},
                  {"v", v_1_out},
                  {"a", a_1_out}}},
                {"coupling_2_outer", {
                  {"id", "coupling_2_outer"},
                  {"gamma", 2},
                  {"x", solver_with_gap.get_points()[2]},
                  {"p", p_2_out},
                  {"q", q_2_out},
                  {"r", r_2_out},
                  {"v", v_2_out},
                  {"a", a_2_out}}},
                {"coupling_2_inner", {
                  {"id", "coupling_2_inner"},
                  {"gamma", 2},
                  {"x", solver_with_gap.get_points()[3]},
                  {"p", p_2_in},
                  {"q", q_2_in},
                  {"r", r_2_in},
                  {"v", v_2_in},
                  {"a", a_2_in}}},
                }
              },
              {"time", t},
              {"units", {
                {"p", "g s^{-2} cm^{-1}"},
                {"q", "cm^{3} s^{-1}"},
                {"a", "cm^{2}"},
                {"r", "cm"},
                {"v", "cm s^{-1}"}}},
              {"vessel_ids", {
                "coupling_1_inner",
                "coupling_1_outer",
                "coupling_2_outer",
                "coupling_2_inner",
              }}
            };

            f << j.dump(1);

          */
        }

      // output every 100
      if ((i+1) % int(1e-1/tau) == 0) {
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
