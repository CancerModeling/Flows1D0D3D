////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner, Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <fmt/format.h>
#include <utility>

#include "petsc.h"

#include "../cmake-build-debug/_deps/json-src/single_include/nlohmann/json.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/heart_to_breast_1d_solver.hpp"
#include "macrocirculation/heart_to_breast_3d_solver.hpp"
#include "macrocirculation/libmesh_utils.hpp"
#include "macrocirculation/quantities_of_interest.hpp"


namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

std::vector< mc::VesselTipCurrentCouplingData > read_coupling_data(const std::string& filepath)
{

  std::vector< mc::VesselTipCurrentCouplingData > coupling_data;

  using json = nlohmann::json;

  std::fstream file(filepath, std::ios::in);
  json j;
  file >> j;

  for (auto& v: j["vertices"])
  {
    mc::VesselTipCurrentCouplingData data { {0,0,0}, 0, 0, 0, 0, 0, 0 };
    data.p = {v["p"][0], v["p"][1], v["p"][2]};
    data.vertex_id = v["vertex_id"];
    data.pressure = v["pressure"];
    data.concentration = v["concentration"];
    data.R2 = v["R2"];
    data.radius = v["radius"];
    data.level = v["level"];

    coupling_data.push_back(data);
  }

  return coupling_data;
}

int main(int argc, char *argv[]) {
  // Libmesh init
  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  char petsc_version[1000];
  PetscGetVersion(petsc_version, 1000);
  std::cout << petsc_version << std::endl;

  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                                                 //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_3d_model/"))         //
    ("input-pressures-file", "file for the input pressures", cxxopts::value<std::string>()->default_value("./data/3d-input-pressures/linearized_1d0d3d.json"))         //
    ("input-file", "input filename for parameters", cxxopts::value<std::string>()->default_value(""))                                            //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.02"))                                                         //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/3d-meshes/test_full_1d0d3d_cm.e"))                //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const auto out_dir = args["output-directory"].as<std::string>();

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  // setup 3D solver
  log("setting up 3D solver\n");
  auto filename = args["input-file"].as<std::string>();
  // read input parameters
  mc::HeartToBreast3DSolverInputDeck input(filename);
  if (filename.empty()) {
    input.d_T = 1.;
    input.d_dt = 0.1;
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
  log("creating model\n");
  auto solver_3d = mc::HeartToBreast3DSolver(MPI_COMM_WORLD, comm,
                                             input, mesh, eq_sys, p_cap, p_tis,
                                             nut_cap, nut_tis, tum,
                                             K_tis, Dnut_tis_field,
                                             N_bar_cap_field, N_bar_surf_cap_field,
                                             log);
  eq_sys.init();
  solver_3d.set_conductivity_fields();

  auto coupling_data = read_coupling_data(args["input-pressures-file"].as<std::string>());
  solver_3d.setup_1d3d(coupling_data);

  // finalize 3D solver setup
  log("finalizing setup of 3D solver\n");
  solver_3d.setup();
  solver_3d.write_output();

  // NOTE to get relevant values from 3D system to solve 1D system
  // call get_vessel_tip_data_3d()
  // data_3d contains vector of coefficients a and b and also weighted avg of 3D pressure
  auto data_3d = solver_3d.get_vessel_tip_data_3d();
  log("initial 3D data at outlet tips");
  for (const auto &a : data_3d)
    log(fmt::format("avg_p = {}, avg_nut = {}\n", a.d_p_3d_w, a.d_nut_3d_w));


  // Some condition to solve the 3D system
  {
    log("update 1d data in 3d solver\n");
    solver_3d.update_1d_data(coupling_data);

    log("solve 3d systems\n");
    solver_3d.solve();

    log("update 3d data for 1d systems\n");
    solver_3d.update_3d_data();

    solver_3d.write_output();

    // recompute avg 3d values at outlet tips
    data_3d = solver_3d.get_vessel_tip_data_3d();
  }
}