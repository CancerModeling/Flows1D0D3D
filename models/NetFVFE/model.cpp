////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "model.hpp"
#include "rw/vtk_io.hpp"

namespace {

steady_clock::time_point clock_begin = steady_clock::now();
steady_clock::time_point clock_end = steady_clock::now();

void reset_clock() {
  clock_begin = steady_clock::now();
}

std::ostringstream oss;

void random_init() { srand(time(nullptr)); }

void log_msg(std::string &msg, util::Logger &log) {
  if (!msg.empty()) {
    log(msg, "debug");
    msg.clear();
  }
}

} // namespace

namespace netfvfe {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(netfvfe::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Tumor")
    sys.project_solution(netfvfe::initial_condition_tum, nullptr, es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(netfvfe::initial_condition_hyp, nullptr, es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(netfvfe::initial_condition_nec, nullptr, es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(netfvfe::initial_condition_taf, nullptr, es.parameters);
  else if (system_name == "ECM")
    sys.project_solution(netfvfe::initial_condition_ecm, nullptr, es.parameters);
  else if (system_name == "MDE")
    sys.project_solution(netfvfe::initial_condition_mde, nullptr, es.parameters);
  else if (system_name == "Pressure")
    sys.project_solution(netfvfe::initial_condition_pres, nullptr,
                         es.parameters);
  else {
    return;
  }
}
} // namespace netfvfe

// Model setup and run
void netfvfe::model_setup_run(int argc, char **argv,
                             const std::string &filename,
                             Parallel::Communicator *comm) {

  reset_clock();

  // init seed for random number
  random_init();

  // read input file
  auto input = InpDeck(filename);

  // create logger
  std::string logger_file = input.d_output_path + "info";
  if (!input.d_outfile_tag.empty())
    logger_file += "_" + input.d_outfile_tag;
  util::Logger log(logger_file, comm, !input.d_quiet);

  // disable reference counter information
  if (input.d_quiet)
    ReferenceCounter::disable_print_counter_info();


  log("Model: NetFVFE\n", "init");

  // create mesh
  oss << "Setup [Mesh] -> ";
  ReplicatedMesh mesh(*comm);
  if (input.d_read_mesh_flag)
    mesh.read(input.d_mesh_filename);
  else
    util::create_mesh(input, mesh);

  oss << " [Tumor system] -> ";
  EquationSystems tum_sys(mesh);
  // add parameters to system
  tum_sys.parameters.set<InpDeck *>("input_deck") = &input;
  tum_sys.parameters.set<Real>("time_step") = input.d_dt;

  // read if available
  if (input.d_restart)
    tum_sys.read(input.d_sol_restart_file, READ);

  // Add systems, variables and assemble
  auto &nut = tum_sys.add_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = tum_sys.add_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = tum_sys.add_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = tum_sys.add_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = tum_sys.add_system<TransientLinearImplicitSystem>("TAF");
  auto &ecm = tum_sys.add_system<TransientLinearImplicitSystem>("ECM");
  auto &mde = tum_sys.add_system<TransientLinearImplicitSystem>("MDE");
  auto &pres = tum_sys.add_system<TransientLinearImplicitSystem>("Pressure");
  auto &grad_taf =
      tum_sys.add_system<TransientLinearImplicitSystem>("TAF_Gradient");
  auto &vel =
      tum_sys.add_system<TransientLinearImplicitSystem>("Velocity");

  // some initial setups
  {
    if (input.d_restart) {
      nut.update();
      tum.update();
      hyp.update();
      nec.update();
      taf.update();
      ecm.update();
      mde.update();
      pres.update();
      grad_taf.update();
      vel.update();
    }

    if (!input.d_restart) {
      // variable in nutrient system
      nut.add_variable("nutrient", CONSTANT, MONOMIAL);

      // variable in tumor system
      tum.add_variable("tumor", FIRST);
      tum.add_variable("chemical_tumor", FIRST);

      // variable in hypoxic system
      hyp.add_variable("hypoxic", FIRST);

      // variable in necrotic system
      nec.add_variable("necrotic", FIRST);

      // variable in TAF system
      taf.add_variable("taf", FIRST);

      // variable in ECM system
      ecm.add_variable("ecm", FIRST);

      // variable in MDE system
      mde.add_variable("mde", FIRST);

      // variable in Pressure system
      pres.add_variable("pressure", CONSTANT, MONOMIAL);

      // variable in gradient TAF system
      grad_taf.add_variable("taf_gradx", FIRST);
      grad_taf.add_variable("taf_grady", FIRST);
      if (input.d_dim > 2)
        grad_taf.add_variable("taf_gradz", FIRST);

      // variable in velocity system
      vel.add_variable("velocity_x", FIRST);
      vel.add_variable("velocity_y", FIRST);
      if (input.d_dim > 2)
        vel.add_variable("velocity_z", FIRST);

      // attach initial condition function to systems
      tum.attach_init_function(initial_condition);
      hyp.attach_init_function(initial_condition);
      nut.attach_init_function(initial_condition);
      ecm.attach_init_function(initial_condition);
      mde.attach_init_function(initial_condition);
      pres.attach_init_function(initial_condition);
    }

    // Add boundary condition
    boundary_condition_nut(tum_sys);
  }

  //
  // Create Model class
  //
  auto model = Model(argc, argv, filename, comm, input, mesh, tum_sys,
                     nec, tum, nut, hyp, taf, ecm, mde, pres, grad_taf, vel,
                     log);

  // run model
  model.run();
}

// Model class
netfvfe::Model::Model(
    int argc, char **argv, const std::string &filename, Parallel::Communicator *comm,
    InpDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
    TransientLinearImplicitSystem &nec, TransientLinearImplicitSystem &tum,
    TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &hyp,
    TransientLinearImplicitSystem &taf, TransientLinearImplicitSystem &ecm,
    TransientLinearImplicitSystem &mde, TransientLinearImplicitSystem &pres,
    TransientLinearImplicitSystem &grad_taf,
    TransientLinearImplicitSystem &vel,
    util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log),
      d_network(this),
      d_nec_assembly(this, "Necrotic", d_mesh, nec),
      d_tum_assembly(this, "Tumor", d_mesh, tum),
      d_nut_assembly(this, "Nutrient", d_mesh, nut),
      d_hyp_assembly(this, "Hypoxic", d_mesh, hyp),
      d_taf_assembly(this, "TAF", d_mesh, taf),
      d_ecm_assembly(this, "ECM", d_mesh, ecm),
      d_mde_assembly(this, "MDE", d_mesh, mde),
      d_pres_assembly(this, "Pressure", d_mesh, pres),
      d_grad_taf_assembly(this, "TAF_Gradient", d_mesh, grad_taf),
      d_vel_assembly(this, "Velocity", d_mesh, vel), d_ghosting_fv(d_mesh) {

  d_nut_id = d_nut_assembly.d_sys.number();
  d_tum_id = d_tum_assembly.d_sys.number();
  d_hyp_id = d_hyp_assembly.d_sys.number();
  d_nec_id = d_nec_assembly.d_sys.number();
  d_taf_id = d_taf_assembly.d_sys.number();
  d_ecm_id = d_ecm_assembly.d_sys.number();
  d_mde_id = d_mde_assembly.d_sys.number();
  d_pres_id = d_pres_assembly.d_sys.number();
  d_grad_taf_id = d_grad_taf_assembly.d_sys.number();
  d_vel_id = d_vel_assembly.d_sys.number();
  d_pres_1d_id = d_vel_id + 1;
  d_nut_1d_id = d_vel_id + 2;

  d_sys_names.resize(d_vel_id + 3);
  d_sys_names[d_nut_id] = "Nutrient";
  d_sys_names[d_tum_id] = "Tummor";
  d_sys_names[d_hyp_id] = "Hypoxic";
  d_sys_names[d_nec_id] = "Necrotic";
  d_sys_names[d_taf_id] = "TAF";
  d_sys_names[d_ecm_id] = "ECM";
  d_sys_names[d_mde_id] = "MDE";
  d_sys_names[d_pres_id] = "Pressure";
  d_sys_names[d_grad_taf_id] = "TAF_Gradient";
  d_sys_names[d_vel_id] = "Velocity";
  d_sys_names[d_pres_1d_id] = "Pressure_1D";
  d_sys_names[d_nut_1d_id] = "Nutrient_1D";

  // init timestep log
  d_log.init_ts(d_sys_names);

  // bounding box
  d_bounding_box.first =
      Point(d_input.d_domain_params[0], d_input.d_domain_params[2],
            d_input.d_domain_params[4]);
  d_bounding_box.second =
      Point(d_input.d_domain_params[1], d_input.d_domain_params[3],
            d_input.d_domain_params[5]);

  // remaining system setup
  {
    // attach assembly objects to various systems
    nec.attach_assemble_object(d_nec_assembly);
    tum.attach_assemble_object(d_tum_assembly);
    nut.attach_assemble_object(d_nut_assembly);
    hyp.attach_assemble_object(d_hyp_assembly);
    taf.attach_assemble_object(d_taf_assembly);
    ecm.attach_assemble_object(d_ecm_assembly);
    mde.attach_assemble_object(d_mde_assembly);
    pres.attach_assemble_object(d_pres_assembly);
    grad_taf.attach_assemble_object(d_grad_taf_assembly);
    vel.attach_assemble_object(d_vel_assembly);

    // add ghosting functors
    //pres.get_dof_map().add_coupling_functor(d_ghosting_fv);
    //nut.get_dof_map().add_coupling_functor(d_ghosting_fv);

    //
    // Initialize and print system
    //
    if (!d_input.d_restart)
      d_tum_sys.init();

    nut.time = d_input.d_init_time;
    tum.time = d_input.d_init_time;
    hyp.time = d_input.d_init_time;
    nec.time = d_input.d_init_time;
    taf.time = d_input.d_init_time;
    ecm.time = d_input.d_init_time;
    mde.time = d_input.d_init_time;
    pres.time = d_input.d_init_time;
    grad_taf.time = d_input.d_init_time;
    vel.time = d_input.d_init_time;

    //    // set Petsc matrix option to suppress the error
    //    {
    //      PetscMatrix<Number> *pet_mat =
    //          dynamic_cast<PetscMatrix<Number> *>(pres.matrix);
    //      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    //    }
    //    {
    //      PetscMatrix<Number> *pet_mat =
    //          dynamic_cast<PetscMatrix<Number> *>(nut.matrix);
    //      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    //    }

    if (d_input.d_perform_output and !d_input.d_quiet) {
      d_delayed_msg += "Libmesh Info \n";
      d_delayed_msg += d_mesh.get_info();
      d_delayed_msg += " \n";
      d_delayed_msg +=d_tum_sys.get_info();
      d_mesh.write("mesh_" + d_input.d_outfile_tag + ".e");
    }
  }

  // 1-D network
  oss << "[Network] \n";
  d_log(oss, "init");
  d_network.create_initial_network();
  d_log(d_delayed_msg, "debug");

  // initialize qoi data
  d_qoi = util::QoIVec({"tumor_mass", "hypoxic_mass", "necrotic_mass",
                        "tumor_l2"});

  // save setup end time
  clock_end = steady_clock::now();
  d_log.d_setup_time = util::TimePair(clock_begin, clock_end);
}

void netfvfe::Model::run() {

  // print initial state
  {
    // Tumor model
    if (d_input.d_perform_output)
      write_system(0);

    // network
    d_network.writeDataToVTKTimeStep_VGM(0);
  }

  // set time parameters
  d_step = d_input.d_init_step;
  d_time = d_input.d_init_time;
  d_dt = d_input.d_dt;

  d_tum_sys.parameters.set<unsigned int>("linear solver maximum iterations") =
      d_input.d_linear_max_iters;

  //
  // Solve step
  //

  if (!d_input.d_test_name.empty()) {
    oss << " \n Solving sub-system: " << d_input.d_test_name << " \n";
    d_log(oss, "debug");
  }

  // check for tumor-network test
  if (d_input.d_test_name == "test_net_tum_2" or
      d_input.d_test_name == "test_nut" or
      d_input.d_test_name == "test_nut_2" or
      d_input.d_test_name == "test_pressure") {

    // solve for pressure only once
    auto nt = d_input.d_nonlin_max_iters;
    d_input.d_nonlin_max_iters = 2 * nt;
    solve_pressure();
    d_input.d_nonlin_max_iters = nt;
  }

  // check for pressure test
  if (d_input.d_test_name == "test_pressure") {

    // write tumor solution
    write_system(1);
    d_network.writeDataToVTKTimeStep_VGM(1);
    return;
  }

  // based on test_name, decide what systems to output
  if (!d_input.d_test_name.empty()) {

    // hide output to mde, ecm, grad taf
    d_mde_assembly.d_sys.hide_output() = true;
    d_ecm_assembly.d_sys.hide_output() = true;

    if (d_input.d_test_name == "test_nut" or
        d_input.d_test_name == "test_nut_2" or
        d_input.d_test_name == "test_pressure") {
      d_taf_assembly.d_sys.hide_output() = true;
      d_grad_taf_assembly.d_sys.hide_output() = true;
    }
  }

  do {

    // Prepare time step
    d_step++;
    d_time += d_dt;

    // init ts log
    d_log.ready_new_step(int(d_step) - 1);
    auto solve_clock = steady_clock::now();

    // check if this is output step
    d_is_output_step = false;
    if ((d_step >= d_input.d_dt_output_interval and
         (d_step % d_input.d_dt_output_interval == 0)) or
        d_step == 0)
      d_is_output_step = true;

    d_is_growth_step = false;
    if (d_step % d_input.d_network_update_interval == 0)
      d_is_growth_step = true;

    oss << "Time step: " << d_step << ", time: " << d_time << "\n";
    d_log(oss, "integrate");
    d_log(" \n", "integrate");

    // solve tumor-network system
    solve_system();

    // compute qoi
    compute_qoi();

    // update network
    if (d_is_growth_step) {
      d_log("  Updating Network\n", "net update");
      d_network.updateNetwork(d_taf_assembly, d_grad_taf_assembly);
      d_log(" \n", "net update");
    }

    // Post-processing
    if (d_is_output_step) {

      // write tumor solution
      write_system((d_step - d_input.d_init_step) /
                       d_input.d_dt_output_interval);
      d_network.writeDataToVTKTimeStep_VGM((d_step - d_input.d_init_step) /
                                           d_input.d_dt_output_interval);
    }

    // output qoi
    if (d_step == 1)
      d_log.log_qoi_header(d_time, d_qoi.get_last(), d_qoi.get_names());
    else
      d_log.log_qoi(d_time, d_qoi.get_last());

    // add to log
    d_log.add_solve_time(util::TimePair(solve_clock, steady_clock::now()));

    // output time logger info
    d_log.log_ts();

  } while (d_step < d_input.d_steps);
}

void netfvfe::Model::write_system(const unsigned int &t_step) {

  //ExodusII_IO exodus(d_mesh);

  // write mesh and simulation results
  //std::string filename = d_input.d_outfilename + ".e";

  // write to exodus
  //exodus.write_timestep(filename, d_tum_sys, 1, d_time);

  //
  rw::VTKIO(d_mesh).write_equation_systems(
      d_input.d_outfilename + "_" + std::to_string(t_step) + ".pvtu",
      d_tum_sys);

  // save for restart
  if (d_input.d_restart_save &&
      (d_step % d_input.d_dt_restart_save_interval == 0)) {

    oss << "\n  Saving files for restart\n";
    d_log(oss);
    if (t_step == 0) {

      std::string mesh_file = "mesh_" + d_input.d_outfile_tag + ".e";
      d_mesh.write(mesh_file);
    }

    std::string solutions_file = "solution_" + d_input.d_outfile_tag + "_" +
                                 std::to_string(d_step) + ".e";
    d_tum_sys.write(solutions_file, WRITE);
  }

  // debug
  // output 3D pressure and nutrient stored in network to txt file
  // header
  // pressure nutrient taf
  //if (d_comm_p->size() > 1)
    util::localize_solution_with_elem_id_numbering_non_const_elem(
      d_taf_assembly, d_network.localized_taf_3D, d_network.phi_TAF, {0},
      false);

  if (d_comm_p->rank() == 0 and d_step > 1) {

    std::ofstream of;
    of.open(d_input.d_outfilename + "_debug_" + std::to_string(t_step) +
            ".txt");
    for (unsigned int i = 0; i < d_network.N_tot_3D; i++) {

      double pres_val = 0.;
      double nut_val = 0.;
      double taf_val = 0.;

      if (d_comm_p->size() == 1) {

        const auto *elem = d_mesh.elem_ptr(i);
        d_pres_assembly.init_dof(elem);
        d_nut_assembly.init_dof(elem);
        d_taf_assembly.init_dof(elem);

        pres_val = d_pres_assembly.get_current_sol(0);
        nut_val = d_nut_assembly.get_current_sol(0);
        taf_val = 0.;
        for (unsigned int j=0; j < d_taf_assembly.d_phi.size(); j++)
          taf_val += d_taf_assembly.get_current_sol(j);
        taf_val = taf_val / d_taf_assembly.d_phi.size();

      } else {

        pres_val = d_network.P_3D[i];
        nut_val = d_network.phi_sigma_3D[i];
        taf_val = d_network.phi_TAF[i];
      }

      of << pres_val << " " << nut_val << " " << taf_val << "\n";
    }
    of.close();
  }
}

void netfvfe::Model::compute_qoi() {

  // integral of total tumor
  double value_mass = 0.;
  double total_mass = 0.;
  util::computeMass(d_tum_sys, "Tumor", "tumor", value_mass);
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  // add to qoi data
  double tumor_mass = total_mass;

  // integral of hypoxic
  value_mass = 0.; total_mass = 0.;
  util::computeMass(d_tum_sys, "Hypoxic", "hypoxic", value_mass);
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  // add to qoi data
  double hypoxic_mass = total_mass;

  // integral of necrotic
  value_mass = 0.; total_mass = 0.;
  util::computeMass(d_tum_sys, "Necrotic", "necrotic", value_mass);
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  // add to qoi data
  double necrotic_mass = total_mass;

  // L2 norm of total tumor
  value_mass = 0.; total_mass = 0.;
  value_mass = d_tum_assembly.d_sys.solution->l2_norm();
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  // add to qoi data
  double tumor_L2_norm = total_mass;

  d_qoi.add({tumor_mass, hypoxic_mass, necrotic_mass, tumor_L2_norm});
}

void netfvfe::Model::solve_system() {

  if (!d_input.d_test_name.empty()) {
    if (d_input.d_test_name == "test_nut")
      test_nut();
    else if (d_input.d_test_name == "test_nut_2")
      test_nut_2();
    else if (d_input.d_test_name == "test_taf")
      test_taf();
    else if (d_input.d_test_name == "test_taf_2")
      test_taf_2();
    else if (d_input.d_test_name == "test_tum")
      test_tum();
    else if (d_input.d_test_name == "test_tum_2")
      test_tum_2();
    else if (d_input.d_test_name == "test_net_tum")
      test_net_tum();
    else if (d_input.d_test_name == "test_net_tum_2")
      test_net_tum_2();
    else
      libmesh_error_msg("Test name " + d_input.d_test_name +
                        " not recognized.");

    return;
  }

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &ecm = d_tum_sys.get_system<TransientLinearImplicitSystem>("ECM");
  auto &mde = d_tum_sys.get_system<TransientLinearImplicitSystem>("MDE");
  auto &pres = d_tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  ecm.time = d_time;
  mde.time = d_time;
  pres.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;
  *taf.old_local_solution = *taf.current_local_solution;
  *ecm.old_local_solution = *ecm.current_local_solution;
  *mde.old_local_solution = *mde.current_local_solution;
  *pres.old_local_solution = *pres.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // reset nonlinear step
  d_nonlinear_step = 0;

  // check if we are decoupling the nutrients
  if (d_input.d_decouple_nutrients) {

    reset_clock();

    d_log("      Solving |1D nutrient|\n", "solve sys");
    d_log( " \n", "solve sys");
    d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

    d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
  }

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solver for 1D nutrient
    if (!d_input.d_decouple_nutrients) {
      reset_clock();

      d_log("|1D nutrient| -> ", "solve sys");
      d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

      d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
    }

    // solve nutrient
    reset_clock();
    d_log("|3D nutrient| -> ", "solve sys");
    nut.solve();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);

    // solve tumor
    reset_clock();
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    d_log("|tumor| -> ", "solve sys");
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();
    d_log.add_sys_solve_time(clock_begin, d_tum_id);

    // solve hypoxic
    reset_clock();
    d_log("|hypoxic| -> ", "solve sys");
    hyp.solve();
    d_log.add_sys_solve_time(clock_begin, d_hyp_id);

    // solve necrotic
    reset_clock();
    d_log("|necrotic| -> ", "solve sys");
    nec.solve();
    d_log.add_sys_solve_time(clock_begin, d_nec_id);

    // solve taf
    reset_clock();
    d_log("|taf| -> ", "solve sys");
    taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_taf_id);

    // solve mde
    reset_clock();
    d_log("|mde| -> ", "solve sys");
    mde.solve();
    d_log.add_sys_solve_time(clock_begin, d_mde_id);

    // solve ecm
    reset_clock();
    d_log("|ecm|\n", "solve sys");
    ecm.solve();
    d_log.add_sys_solve_time(clock_begin, d_ecm_id);

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);

  // solve for gradient of taf
  reset_clock();
  d_log("      Solving |grad taf| \n", "solve sys");
  grad_taf.solve();
  d_log.add_sys_solve_time(clock_begin, d_grad_taf_id);
}

void netfvfe::Model::solve_pressure() {

  auto solve_clock = steady_clock::now();
  reset_clock();

  // get systems
  auto &pres = d_tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");
  auto &vel = d_tum_sys.get_system<TransientLinearImplicitSystem>("Velocity");

  // update time
  pres.time = d_time;
  vel.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *pres.old_local_solution = *pres.current_local_solution;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_pres(
      pres.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // reset nonlinear step
  d_nonlinear_step = 0;

  // Solve pressure system
  d_log("  Nonlinear loop for pressure\n\n", "solve pres");
  d_log(" \n", "solve pres");
  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
    // Debug
    // for (unsigned int l = 0; l < 10; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l;
    d_log(oss, "solve pres");

    // solver for 1-D pressure and nutrient
    reset_clock();
    d_log( " |1D pressure| -> ", "solve pres");
    d_network.solveVGMforPressure(d_pres_assembly);
    if (d_log.d_cur_step >= 0)
      d_log.add_sys_solve_time(clock_begin, d_pres_1d_id);


    // solve for pressure in tissue
    reset_clock();
    last_nonlinear_soln_pres->zero();
    last_nonlinear_soln_pres->add(*pres.solution);
    d_log("|3D pressure|\n", "solve pres");
    pres.solve();

    last_nonlinear_soln_pres->add(-1., *pres.solution);
    last_nonlinear_soln_pres->close();
    if (d_log.d_cur_step >= 0)
      d_log.add_sys_solve_time(clock_begin, d_pres_id);

    // Nonlinear iteration error
    double nonlinear_loc_error_pres = last_nonlinear_soln_pres->linfty_norm();
    double nonlinear_global_error_pres = 0.;
    MPI_Allreduce(&nonlinear_loc_error_pres, &nonlinear_global_error_pres, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = pres.n_linear_iterations();
      const Real final_linear_residual = pres.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error_pres << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error_pres < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "\n      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge pres");

      break;
    }

    log_msg(d_delayed_msg, d_log);

  } // nonlinear solver loop

  d_log(" \n", "solve pres");

  clock_end = std::chrono::steady_clock::now();
  if (d_log.d_cur_step >= 0) {
    d_log.add_pres_solve_time(util::TimePair(solve_clock, clock_end));
    d_log.add_pres_nonlin_iter(d_nonlinear_step);
  }

  // solve for velocity
  reset_clock();
  d_log("      Solving |velocity| \n", "solve pres");
  d_log(" \n", "solve pres");
  vel.solve();
  if (d_log.d_cur_step >= 0)
    d_log.add_sys_solve_time(clock_begin, d_vel_id);
}

void netfvfe::Model::test_nut() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

  // update time
  nut.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // reset nonlinear step
  d_nonlinear_step = 0;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_nut(
      nut.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");


  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
    // for (unsigned int l = 0; l < 1; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solver for 1-D nutrient
    reset_clock();
    d_log("|1D nutrient| -> ", "solve sys");
    d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);
    d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);

    // solve for nutrient in tissue
    last_nonlinear_soln_nut->zero();
    last_nonlinear_soln_nut->add(*nut.solution);
    reset_clock();
    d_log("|3D nutrient|\n", "solve sys");
    nut.solve();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);
    last_nonlinear_soln_nut->add(-1., *nut.solution);
    last_nonlinear_soln_nut->close();

    double nonlinear_loc_error_nut = last_nonlinear_soln_nut->linfty_norm();
    double nonlinear_global_error_nut = 0.;
    MPI_Allreduce(&nonlinear_loc_error_nut, &nonlinear_global_error_nut, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      // const unsigned int n_linear_iterations = tum.n_linear_iterations();
      // const Real final_linear_residual = tum.final_linear_residual();
      const unsigned int n_linear_iterations = nut.n_linear_iterations();
      const Real final_linear_residual = nut.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error_nut << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error_nut < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }

  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);
}

void netfvfe::Model::test_nut_2() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

  // update time
  nut.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // solver for 1-D nutrient
  reset_clock();
  d_log("      Solving |1D nutrient| -> ", "solve sys");
  d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);
  d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);

  reset_clock();
  d_log("|3D nutrient| \n", "solve sys");
  nut.solve();
  d_log.add_sys_solve_time(clock_begin, d_nut_id);
}

void netfvfe::Model::test_taf() {

  // get systems
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  taf.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *taf.old_local_solution = *taf.current_local_solution;

  reset_clock();
  d_log("      Solving |taf| \n", "solve sys");
  taf.solve();
  d_log.add_sys_solve_time(clock_begin, d_taf_id);

  // solve for gradient of taf
  reset_clock();
  d_log("|grad taf| \n", "solve sys");
  grad_taf.solve();
  d_log.add_sys_solve_time(clock_begin, d_grad_taf_id);
}

void netfvfe::Model::test_taf_2() {

  // get systems
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  taf.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // solve for pressure
  solve_pressure();

  // update old solution
  *taf.old_local_solution = *taf.current_local_solution;

  reset_clock();
  d_log("      Solving |taf| \n", "solve sys");
  taf.solve();
  d_log.add_sys_solve_time(clock_begin, d_taf_id);

  // solve for gradient of taf
  reset_clock();
  d_log("|grad taf| \n", "solve sys");
  grad_taf.solve();
  d_log.add_sys_solve_time(clock_begin, d_grad_taf_id);
}

void netfvfe::Model::test_tum() {

  // set nutrient value assuming cylindrical source
  {
    double L = d_input.d_domain_params[1];
    double R = 0.05 * L;
    Point xc = Point(L - 2. * R, L - 2. * R, 0.);

    // Looping through elements
    for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

      d_nut_assembly.init_dof(elem);

      auto x = elem->centroid();
      Point dx = Point(x(0), x(1), 0.) - xc;
      if (dx.norm() < R) {
        auto nut = d_nut_assembly.get_current_sol(0);
        if (nut < 1.)
          d_nut_assembly.d_sys.solution->set(
              d_nut_assembly.get_global_dof_id(0), 1.);
      }
    }

    d_nut_assembly.d_sys.solution->close();
    d_nut_assembly.d_sys.update();
  }

  // get systems
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");

  // update time
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;

  // reset nonlinear step
  d_nonlinear_step = 0;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solve tumor
    reset_clock();
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    d_log("|tumor| -> ", "solve sys");
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();
    d_log.add_sys_solve_time(clock_begin, d_tum_id);

    // solve hypoxic
    reset_clock();
    d_log("|hypoxic| -> ", "solve sys");
    hyp.solve();
    d_log.add_sys_solve_time(clock_begin, d_hyp_id);

    // solve necrotic
    reset_clock();
    d_log("|necrotic|\n", "solve sys");
    nec.solve();
    d_log.add_sys_solve_time(clock_begin, d_nec_id);

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);
}

void netfvfe::Model::test_tum_2() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;

  // reset nonlinear step
  d_nonlinear_step = 0;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solve nutrient
    reset_clock();
    d_log("|3D nutrient| -> ", "solve sys");
    nut.solve();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);

    // solve tumor
    reset_clock();
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    d_log("|tumor| -> ", "solve sys");
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();
    d_log.add_sys_solve_time(clock_begin, d_tum_id);

    // solve hypoxic
    reset_clock();
    d_log("|hypoxic| -> ", "solve sys");
    hyp.solve();
    d_log.add_sys_solve_time(clock_begin, d_hyp_id);

    // solve necrotic
    reset_clock();
    d_log("|necrotic|\n", "solve sys");
    nec.solve();
    d_log.add_sys_solve_time(clock_begin, d_nec_id);

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);
}

void netfvfe::Model::test_net_tum() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // reset nonlinear step
  d_nonlinear_step = 0;

  // check if we are decoupling the nutrients
  if (d_input.d_decouple_nutrients) {
    reset_clock();


    d_log("      Solving |1D nutrient| \n", "solve sys");
    d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

    d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
  }

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solver for 1D nutrient
    if (!d_input.d_decouple_nutrients) {
      reset_clock();

      d_log("|1D nutrient| -> ", "solve sys");
      d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

      d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
    }

    // solve nutrient
    reset_clock();
    d_log("|3D nutrient| -> ", "solve sys");
    nut.solve();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);

    // solve tumor
    reset_clock();
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    d_log("|tumor| -> ", "solve sys");
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();
    d_log.add_sys_solve_time(clock_begin, d_tum_id);

    // solve hypoxic
    reset_clock();
    d_log("|hypoxic| -> ", "solve sys");
    hyp.solve();
    d_log.add_sys_solve_time(clock_begin, d_hyp_id);

    // solve necrotic
    reset_clock();
    d_log("|necrotic|\n", "solve sys");
    nec.solve();
    d_log.add_sys_solve_time(clock_begin, d_nec_id);

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);

  // compute below only when we are performing output as these do not play
  // direct role in evolution of sub-system

  if (d_is_output_step or d_is_growth_step) {

    // solve for taf
    reset_clock();
    *taf.old_local_solution = *taf.current_local_solution;
    d_log("      Solving |taf| -> ", "solve sys");
    taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_taf_id);

    // solve for gradient of taf
    reset_clock();
    d_log("|grad taf|\n", "solve sys");
    grad_taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_grad_taf_id);
  }
  d_log(" \n", "solve sys");
}

void netfvfe::Model::test_net_tum_2() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // reset nonlinear step
  d_nonlinear_step = 0;

  // check if we are decoupling the nutrients
  if (d_input.d_decouple_nutrients) {
    reset_clock();


    d_log("      Solving |1D nutrient| \n", "solve sys");
    d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

    d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
  }

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");
  d_log(" \n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solver for 1D nutrient
    if (!d_input.d_decouple_nutrients) {
      reset_clock();

      d_log("|1D nutrient| -> ", "solve sys");
      d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);

      d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);
    }

    // solve nutrient
    reset_clock();
    d_log("|3D nutrient| -> ", "solve sys");
    nut.solve();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);

    // solve tumor
    reset_clock();
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    d_log("|tumor| -> ", "solve sys");
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();
    d_log.add_sys_solve_time(clock_begin, d_tum_id);

    // solve hypoxic
    reset_clock();
    d_log("|hypoxic| -> ", "solve sys");
    hyp.solve();
    d_log.add_sys_solve_time(clock_begin, d_hyp_id);

    // solve necrotic
    reset_clock();
    d_log("|necrotic|\n", "solve sys");
    nec.solve();
    d_log.add_sys_solve_time(clock_begin, d_nec_id);

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);

  // compute below only when we are performing output as these do not play
  // direct role in evolution of sub-system

  if (d_is_output_step or d_is_growth_step) {

    // solve for taf
    reset_clock();
    *taf.old_local_solution = *taf.current_local_solution;
    d_log("      Solving |taf| -> ", "solve sys");
    taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_taf_id);

    // solve for gradient of taf
    reset_clock();
    d_log("|grad taf|\n", "solve sys");
    grad_taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_grad_taf_id);
  }
  d_log(" \n", "solve sys");
}