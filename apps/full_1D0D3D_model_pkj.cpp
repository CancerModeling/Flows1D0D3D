////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner, Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <fmt/format.h>
#include <macrocirculation/graph_csv_writer.hpp>
#include <memory>
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
#include "macrocirculation/heart_to_breast_3d_solver.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/libmesh_utils.hpp"
#include "macrocirculation/nonlinear_linear_coupling.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {

  // Libmesh init
  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                                         //
    ("tau-1d", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))           //
    ("tau-3d", "time step size for the 3D model", cxxopts::value<double>()->default_value("1."))                               //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("100."))                             //
    ("tau-out", "Simulation output interval", cxxopts::value<double>()->default_value("5."))                             //
    ("t-3d-start", "Simulation start time for 3D model", cxxopts::value<double>()->default_value("10."))                             //
    ("n-1d-solves", "Number of times 1D equation is solved per macroscale time", cxxopts::value<int>()->default_value("10"))                   //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_full_1d0d3d_pkj/")) //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.02"))                                                 //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/meshes/test_full_1d0d3d_cm.e"))           //
    ("deactivate-3d-1d-coupling", "deactivates the 3d-1d coupling", cxxopts::value<bool>()->default_value("false"))             //
    ("input-file", "input filename for parameters", cxxopts::value<std::string>()->default_value(""))                           //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    // setup 1D solver
    const double t_end = args["t-end"].as<double>();
    const auto t_3d_start = args["t-3d-start"].as<double>();
    const auto n_1d_solves = args["n-1d-solves"].as<int>();
    const auto tau_1d = args["tau-1d"].as<double>();
    const auto tau_3d = args["tau-3d"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();

    const std::size_t max_iter = 160000000;
    auto out_dir = args["output-directory"].as<std::string>();

    const auto activate_3d_1d_coupling = !args["deactivate-3d-1d-coupling"].as<bool>();

    mc::HeartToBreast1DSolver solver_1d(MPI_COMM_WORLD);
    solver_1d.set_output_folder(out_dir);
    solver_1d.setup(degree, tau_1d, mc::BoundaryModel::DiscreteRCRTree);

    // create logger
    mc::Logger log(out_dir + "sim", comm->rank());

    // setup 3D solver
    log("setting up 3D solver\n");
    auto filename = args["input-file"].as<std::string>();
    // read input parameters
    mc::HeartToBreast3DSolverInputDeck input(filename);
    if (filename.empty()) {
      input.d_T = t_end;
      input.d_dt = tau_3d;
      input.d_h = args["mesh-size"].as<double>();
      input.d_mesh_file = args["mesh-file"].as<std::string>();
      input.d_out_dir = out_dir;
      input.d_debug_lvl = 1;
      input.d_perf_regularized = false;
      input.d_perf_fn_type = "const";
      input.d_perf_regularized = true;
      input.d_perf_fn_type = "linear";
      input.d_perf_neigh_size = std::make_pair(4., 10.);
    }
    log("input data \n" + input.print_str() + "\n");

    // create mesh
    log("creating mesh\n");
    lm::ReplicatedMesh mesh(*comm);
    if (!input.d_mesh_file.empty()) {
      mesh.read(input.d_mesh_file);
      //input.d_h = mc::get_min_nodal_spacing(mesh);
      input.d_h = mc::get_mesh_size_estimate_using_element_volume(mesh);
      log(fmt::format("mesh size = {}\n", input.d_h));
    } else {
      long N = long(1. / input.d_h);
      lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                            1., 0., 1., lm::HEX8);
    }

    // create equation system
    log("creating equation system\n");
    lm::EquationSystems eq_sys(mesh);
    eq_sys.parameters.set<mc::HeartToBreast3DSolverInputDeck *>("input_deck") = &input;
    eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
    auto &p_cap = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Capillary_Pressure");
    p_cap.add_variable("p_cap", lm::FIRST);
    auto &p_tis = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tissue_Pressure");
    p_tis.add_variable("p_tis", lm::FIRST);
    auto &nut_cap = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Capillary_Nutrient");
    nut_cap.add_variable("nut_cap", lm::FIRST);
    auto &nut_tis = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tissue_Nutrient");
    nut_tis.add_variable("nut_tis", lm::FIRST);

    // create spatial field of hydraulic conductivity
    auto &K_tis = eq_sys.add_system<lm::ExplicitSystem>("Tissue_K");
    K_tis.add_variable("k_tis", lm::CONSTANT, lm::MONOMIAL);
    auto &Dnut_tis_field = eq_sys.add_system<lm::ExplicitSystem>("Tissue_D_Nut");
    Dnut_tis_field.add_variable("Dnut_tis", lm::CONSTANT, lm::MONOMIAL);
    auto &N_bar_cap_field = eq_sys.add_system<lm::ExplicitSystem>("Avg_Capillary_Surf_Area");
    N_bar_cap_field.add_variable("n_bar_cap", lm::CONSTANT, lm::MONOMIAL);
    auto &N_bar_surf_cap_field = eq_sys.add_system<lm::ExplicitSystem>("Avg_Capillary_Cross_Section_Area");
    N_bar_surf_cap_field.add_variable("n_bar_surf_cap", lm::CONSTANT, lm::MONOMIAL);


    // create model that holds all essential variables
    log("creating model\n");
    auto solver_3d = mc::HeartToBreast3DSolver(MPI_COMM_WORLD, comm,
                       input, mesh, eq_sys, p_cap, p_tis,
                       nut_cap, nut_tis,
                       K_tis, Dnut_tis_field,
                       N_bar_cap_field, N_bar_surf_cap_field,
                       log);
    eq_sys.init();
    solver_3d.set_conductivity_fields();

    // setup the 1D pressure data in 3D solver
    log("setting 1D-3D coupling data in 3D solver\n");
    auto data_1d = solver_1d.get_vessel_tip_pressures();
    solver_3d.setup_1d3d(data_1d);

    // finalize 3D solver setup
    log("finalizing setup of 3D solver\n");
    solver_3d.setup();
    solver_3d.write_output();

    // NOTE to get relevant values from 3D system to solve 1D system
    // call get_vessel_tip_data_3d()
    // data_3d contains vector of coefficients a and b and also weighted avg of 3D pressure
    auto data_3d = solver_3d.get_vessel_tip_data_3d();

    // time integration

    // step 1 - only solve 1d system
    log("jump starting 1d systems\n");
    const auto begin_t = std::chrono::steady_clock::now();
    double t_now = 0.;
    while (t_now <= t_3d_start) {
      solver_1d.solve();
      t_now = solver_1d.get_time();
    }

    // step 2 - solve 3d and 1d systems
    log("starting to solve 1d-3d coupled system\n");
    size_t time_step_3d = 0;
    while (t_now <= t_end) {

      // solve 1d systems n number of times and compute average values at tips
      {
        log("solve 1d system and get average tip data\n");
        data_1d = solver_1d.get_vessel_tip_pressures();
        data_3d = solver_3d.get_vessel_tip_data_3d();

        // set macroscale time as current time in 1d solver
        int n_half = int(n_1d_solves / 2);
        double t_macro = t_now - n_half * tau_1d;
        solver_1d.set_time(t_macro); // FIXME

        // update the boundary conditions of the 1D system:
        log("update 3d average tip values in 1d system\n");
        {
          std::map<size_t, double> new_tip_pressures;

          for (size_t k = 0; k < data_1d.size(); k += 1) {
            auto &d = data_1d[k];
            if (activate_3d_1d_coupling) {
              new_tip_pressures[d.vertex_id] = data_3d.at(k).d_p_3d_w;
            } else {
              // constant 30 mmHg pressures
              new_tip_pressures[d.vertex_id] = 30 * 1.3333;
            }
          }
          solver_1d.update_vessel_tip_pressures(new_tip_pressures);
        }

        // run 1d solver
        log("solve 1d system\n");
        for (int i = -n_half; i < n_half; i++) {
          solver_1d.solve();
          auto data_1d_now = solver_1d.get_vessel_tip_pressures();
          for (size_t j = 0; j < data_1d.size(); j++)
            data_1d[j].pressure += data_1d_now[j].pressure;
        }
        for (auto & j : data_1d)
          j.pressure = j.pressure / (2 * n_half);
      }

      // solve 3d system
      {
        log("update 1d data in 3d solver\n");
        solver_3d.update_1d_data(data_1d);

        log("solve 3d systems\n");
        solver_3d.d_time = t_now;
        solver_3d.solve();

        log("update 3d data for 1d systems\n");
        solver_3d.update_3d_data();
      }

      if (time_step_3d % (int(tau_out / tau_3d)) == 0) {
        solver_3d.write_output();
        solver_1d.write_output();
      }

      // increment time
      t_now += tau_3d;
      time_step_3d++;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}
