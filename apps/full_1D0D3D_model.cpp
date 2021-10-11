
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
#include "macrocirculation/vessel_tree_flow_integrator.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  // Libmesh init
  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  char petsc_version[1000];
  PetscGetVersion(petsc_version, 1000);
  std::cout << petsc_version << std::endl;

  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                                         //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value("1.5625e-05"))                           //
    ("tau-transport", "time step size for the 1D transport", cxxopts::value<double>()->default_value("1.5625e-04"))             //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("5e-1"))                               //
    ("tau-coup", "time step size for updating the coupling", cxxopts::value<double>()->default_value("5e-3"))                   //
    ("t-coup-start", "The time when the 3D coupling gets activated", cxxopts::value<double>()->default_value("2"))              //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("50."))                               //
    ("t-start-averaging-flow", "When should the averaging of the flow start?", cxxopts::value<double>()->default_value("10."))                               //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_full_1d0d3d_pkj/")) //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                                            //
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

  // setup 1D solver
  const double t_end = args["t-end"].as<double>();
  const double t_coup_start = args["t-coup-start"].as<double>();
  const double t_start_averaging_flow = args["t-start-averaging-flow"].as<double>();
  const std::size_t max_iter = 160000000;

  const auto tau = args["tau"].as<double>();
  const auto tau_out = args["tau-out"].as<double>();
  const auto tau_coup = args["tau-coup"].as<double>();
  const auto tau_transport = args["tau-transport"].as<double>();
  auto out_dir = args["output-directory"].as<std::string>();

  if (tau > tau_coup || tau > tau_transport) {
    std::cerr << "flow time step width too large" << std::endl;
    exit(0);
  }

  const auto activate_3d_1d_coupling = !args["deactivate-3d-1d-coupling"].as<bool>();

  const auto output_interval = static_cast<std::size_t>(std::round(tau_out / tau));
  const auto coupling_interval = static_cast<std::size_t>(std::round(tau_coup / tau));
  const auto transport_update_interval = static_cast<std::size_t>(std::round(tau_transport / tau));

  std::cout << output_interval << " " << coupling_interval << " " << transport_update_interval << std::endl;

  if (std::abs(coupling_interval * tau - tau_coup) > 1e-3 * tau) {
    std::cerr << "coupling time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  if (std::abs(transport_update_interval * tau - tau_transport) > 1e-3 * tau) {
    std::cerr << "transport time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  if (std::abs(output_interval * tau - tau_out) > 1e-3 * tau) {
    std::cerr << "output time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  mc::HeartToBreast1DSolver solver_1d(MPI_COMM_WORLD);
  solver_1d.set_output_folder(out_dir);
  solver_1d.setup(degree, tau, mc::BoundaryModel::DiscreteRCRTree);

  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "updating transport every " << transport_update_interval << " interval" << std::endl;

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  // setup 3D solver
  log("setting up 3D solver\n");
  auto filename = args["input-file"].as<std::string>();
  // read input parameters
  mc::HeartToBreast3DSolverInputDeck input(filename);
  if (filename.empty()) {
    input.d_T = t_end;
    input.d_dt = args["time-step"].as<double>();
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
  auto &tum = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tumor");
  tum.add_variable("tum", lm::FIRST);
  tum.add_variable("mu_tum", lm::FIRST);

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
//  log("creating model\n");
//  auto solver_3d = mc::HeartToBreast3DSolver(MPI_COMM_WORLD, comm,
//                                             input, mesh, eq_sys, p_cap, p_tis,
//                                             nut_cap, nut_tis, tum,
//                                             K_tis, Dnut_tis_field,
//                                             N_bar_cap_field, N_bar_surf_cap_field,
//                                             log);
//  eq_sys.init();
//  solver_3d.set_conductivity_fields();

//  // setup the 1D pressure data in 3D solver
//  log("setting 1D-3D coupling data in 3D solver\n");
//  {
//    auto data_1d = solver_1d.get_vessel_tip_pressures();
//    solver_3d.setup_1d3d(data_1d);
//  }
//
//  // finalize 3D solver setup
//  log("finalizing setup of 3D solver\n");
//  solver_3d.setup();
//  solver_3d.write_output();
//
//  // NOTE to get relevant values from 3D system to solve 1D system
//  // call get_vessel_tip_data_3d()
//  // data_3d contains vector of coefficients a and b and also weighted avg of 3D pressure
//  auto data_3d = solver_3d.get_vessel_tip_data_3d();
//  log("initial 3D data at outlet tips");
//  for (const auto &a : data_3d)
//    log(fmt::format("avg_p = {}, avg_nut = {}\n", a.d_p_3d_w, a.d_nut_3d_w));

  double t = 0;
  double t_transport = 0;

  // time integration
  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    // solve flow:
    solver_1d.solve_flow(tau, t);

    if (t-0.5*tau < t_start_averaging_flow && t + 1.5 * tau > t_start_averaging_flow)
      solver_1d.start_0d_flow_integrator();

    t += tau;

    // solve transport:
    if (it % transport_update_interval == 0) {
      solver_1d.solve_transport(tau_transport, t_transport + tau_transport);
      t_transport += tau_transport;
    }

    // solve 3D system:
    if ((t >= t_coup_start - 1e-12) && (it % coupling_interval == 0)) {
      std::cout << "calculates coupling " << std::endl;
      auto data_1d = solver_1d.get_vessel_tip_pressures();

      for (auto &d : data_1d) {
        // just return the values for now:
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "vertex with id = " << d.vertex_id << ", "
                    << "coordinates = (" << d.p.x << ", " << d.p.y << ", " << d.p.z << "), "
                    << "p = " << d.pressure << ", "
                    << "c = " << d.concentration << ", "
                    << "R = " << d.R2 << ", "
                    << "r = " << d.radius << std::endl;
      }

      // Some condition to solve the 3D system
//      {
//        log("update 1d data in 3d solver\n");
//        solver_3d.update_1d_data(data_1d);
//
//        log("solve 3d systems\n");
//        solver_3d.solve();
//
//        log("update 3d data for 1d systems\n");
//        solver_3d.update_3d_data();
//
//        if (it % (20 * coupling_interval) == 0)
//          solver_3d.write_output();
//
//        // recompute avg 3d values at outlet tips
//        data_3d = solver_3d.get_vessel_tip_data_3d();
//      }

      // update the boundary conditions of the 1D system:
//      {
//        std::map<size_t, double> new_tip_pressures;
//        std::map<size_t, double> new_tip_concentrations;
//
//        if (activate_3d_1d_coupling) {
//          std::cout << "size of 3D coupling data is " << data_3d.size() << ", size of 1D coupling data is " << data_1d.size() << std::endl;
//        }
//
//        for (size_t k = 0; k < data_1d.size(); k += 1) {
//          auto &d = data_1d[k];
//          if (activate_3d_1d_coupling) {
//            new_tip_pressures[d.vertex_id] = data_3d.at(k).d_p_3d_w;
//            new_tip_concentrations[d.vertex_id] = data_3d.at(k).d_nut_3d_w; // FIXME: Gets units right
//          } else {
//            // constant 30 mmHg pressures
//            new_tip_pressures[d.vertex_id] = 30 * 1333.3;
//            new_tip_concentrations[d.vertex_id] = 0.;
//          }
//        }
//        solver_1d.update_vessel_tip_pressures(new_tip_pressures);
//        solver_1d.update_vessel_tip_concentrations(new_tip_concentrations);
//      }
    }

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << t << std::endl;

      solver_1d.write_output(t);
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  if (t > t_start_averaging_flow)
    solver_1d.stop_0d_flow_integrator_and_write();

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
}