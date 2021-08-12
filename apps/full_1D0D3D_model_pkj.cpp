////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner, Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <macrocirculation/graph_csv_writer.hpp>
#include <utility>
#include <memory>
#include <fmt/format.h>

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

#include "macrocirculation/heart_to_breast_3d_solver.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {

  // Libmesh init
  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                               //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.))) //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                     //
    ("tau-coup", "time step size for updating the coupling", cxxopts::value<double>()->default_value("1e-3"))         //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("10"))                      //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_full_1d0d3d_pkj/")) //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                                                   //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.02"))                                                         //
    ("hyd-cond-cap", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                                              //
    ("hyd-cond-tis", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                                              //
    ("permeability-cap", "permeability for mass exchange", cxxopts::value<double>()->default_value("0.1"))                                 //
    ("permeability-tis", "permeability for mass exchange", cxxopts::value<double>()->default_value("0.1"))                                 //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value(""))                                                   //
    ("num-points", "number of perfusion points", cxxopts::value<int>()->default_value("10"))                                           //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    // setup 1D solver
    const double t_end = args["t-end"].as<double>();
    const std::size_t max_iter = 160000000;

    const auto tau = args["tau"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();
    const auto tau_coup = args["tau-coup"].as<double>();
    auto out_dir = args["output-directory"].as<std::string>();

    // const double tau_out = tau;
    const auto output_interval = static_cast<std::size_t>(tau_out / tau);
    const auto coupling_interval = static_cast<std::size_t>(tau_coup / tau);

    mc::HeartToBreast1DSolver solver_1d(MPI_COMM_WORLD);
    solver_1d.set_output_folder(out_dir);
    solver_1d.setup(degree, tau, mc::BoundaryModel::DiscreteRCRTree);

    // setup 3D solver
    auto filename = args["input-file"].as<std::string>();
    int num_perf_points = args["num-points"].as<int>();
    // read input parameters
    auto input = mc::HeartToBreast3DSolverInputDeck(filename);
    if (filename == "") {
      input.d_T = t_end;
      input.d_dt = args["time-step"].as<double>();
      input.d_h = args["mesh-size"].as<double>();
      input.d_K_cap = args["hyd-cond-cap"].as<double>();
      input.d_K_tis = args["hyd-cond-tis"].as<double>();
      input.d_Lp_cap = args["permeability-cap"].as<double>();
      input.d_Lp_tis = args["permeability-tis"].as<double>();
      input.d_mesh_file = args["mesh-file"].as<std::string>();
      input.d_out_dir = out_dir;
    }

    // create logger
    mc::Logger log(out_dir + "sim", comm->rank());
    log("input data \n" + input.print_str() + "\n");

    // create mesh
    log("creating mesh\n");
    lm::ReplicatedMesh mesh(*comm);
    long N = long(1. / input.d_h);
    if (input.d_mesh_file != "")
      mesh.read(input.d_mesh_file);
    else
      lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                            1., 0., 1., lm::HEX8);

    // create equation system
    log("creating equation system\n");
    lm::EquationSystems eq_sys(mesh);
    eq_sys.parameters.set<mc::HeartToBreast3DSolverInputDeck *>("input_deck") = &input;
    eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
    auto &p_cap = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Capillary_Pressure");
    p_cap.add_variable("p_cap", lm::FIRST);
    auto &p_tis = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tissue_Pressure");
    p_tis.add_variable("p_tis", lm::FIRST);

    // create spatial field of hydraulic conductivity
    auto &K_cap = eq_sys.add_system<lm::ExplicitSystem>("Capillary_K");
    K_cap.add_variable("k_cap", lm::CONSTANT, lm::MONOMIAL);
    auto &Lp_cap = eq_sys.add_system<lm::ExplicitSystem>("Capillary_Lp");
    Lp_cap.add_variable("lp_cap", lm::CONSTANT, lm::MONOMIAL);

    // create model that holds all essential variables
    log("creating model\n");
    auto solver_3d = mc::HeartToBreast3DSolver(MPI_COMM_WORLD, comm, input, mesh, eq_sys, p_cap, p_tis, K_cap, Lp_cap, log);
    solver_3d.d_dt = input.d_dt;
    solver_3d.setup();

    // time integration
    const auto begin_t = std::chrono::steady_clock::now();
    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver_1d.solve();

      if (it % coupling_interval == 0) {
        std::cout << "calculates coupling " << std::endl;
        auto data = solver_1d.get_vessel_tip_pressures();

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
          solver_1d.update_vessel_tip_pressures(new_tip_pressures);
        }
      }

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", t = " << solver_1d.get_time() << std::endl;

        solver_1d.write_output();
      }

      // break
      if (solver_1d.get_time() > t_end + 1e-12)
        break;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}
