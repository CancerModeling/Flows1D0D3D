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
steady_clock::time_point solve_clock = steady_clock::now();

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

namespace netfcfvfe {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(netfcfvfe::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Prolific")
    sys.project_solution(netfcfvfe::initial_condition_pro, nullptr,
                         es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(netfcfvfe::initial_condition_hyp, nullptr, es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(netfcfvfe::initial_condition_nec, nullptr, es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(netfcfvfe::initial_condition_taf, nullptr, es.parameters);
  else if (system_name == "ECM")
    sys.project_solution(netfcfvfe::initial_condition_ecm, nullptr, es.parameters);
  else if (system_name == "MDE")
    sys.project_solution(netfcfvfe::initial_condition_mde, nullptr, es.parameters);
  else if (system_name == "Pressure")
    sys.project_solution(netfcfvfe::initial_condition_pres, nullptr,
                         es.parameters);
  else {
    return;
  }
}
} // namespace netfcfvfe

// Model setup and run
void netfcfvfe::model_setup_run(int argc, char **argv,
                                const std::string &filename,
                                Parallel::Communicator *comm) {

  auto sim_begin = steady_clock::now();

  reset_clock();

  // init seed for random number
  //random_init();

  // read input file
  auto input = InpDeck(filename);

  input.d_coupled_1d3d = true;

  // create logger
  std::string logger_file = input.d_output_path + "info";
  if (!input.d_outfile_tag.empty())
    logger_file += "_" + input.d_outfile_tag;
  util::Logger log(logger_file, comm, !input.d_quiet);

  // disable reference counter information
  if (input.d_quiet)
    ReferenceCounter::disable_print_counter_info();


  log("Model: NetFCFVFE\n", "init");

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
  auto &pro = tum_sys.add_system<TransientLinearImplicitSystem>("Prolific");
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
  auto &tum =
    tum_sys.add_system<TransientLinearImplicitSystem>("Tumor");

  // some initial setups
  {
    if (input.d_restart) {
      nut.update();
      pro.update();
      hyp.update();
      nec.update();
      taf.update();
      ecm.update();
      mde.update();
      pres.update();
      grad_taf.update();
      vel.update();
      tum.update();
    }

    if (!input.d_restart) {
      // variable in tumor system
      tum.add_variable("tumor", FIRST);
      tum.add_variable("chemical_tumor", FIRST);

      // variable in nutrient system
      nut.add_variable("nutrient", CONSTANT, MONOMIAL);

      // variable in tumor system
      pro.add_variable("prolific", FIRST);
      pro.add_variable("chemical_prolific", FIRST);

      // variable in hypoxic system
      hyp.add_variable("hypoxic", FIRST);
      hyp.add_variable("chemical_hypoxic", FIRST);

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
      pro.attach_init_function(initial_condition);
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
                     nec, pro, nut, hyp, taf, ecm, mde, pres, grad_taf,
                     vel, tum, log);

  // run model
  model.run();

  model.d_log.log_ts_final(
    util::time_diff(sim_begin, steady_clock::now()));
}

// Model class
netfcfvfe::Model::Model(
  int argc, char **argv, const std::string &filename, Parallel::Communicator *comm,
  InpDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
  TransientLinearImplicitSystem &nec, TransientLinearImplicitSystem &pro,
  TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &hyp,
  TransientLinearImplicitSystem &taf, TransientLinearImplicitSystem &ecm,
  TransientLinearImplicitSystem &mde, TransientLinearImplicitSystem &pres,
  TransientLinearImplicitSystem &grad_taf,
  TransientLinearImplicitSystem &vel, TransientLinearImplicitSystem &tum,
  util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log, "NetFCFVFE"),
      d_network(this),
      d_networkVtkWriter(comm, input.d_outfilename_net + "new_"),
      d_nec(this, "Necrotic", d_mesh, nec),
      d_pro(this, "Prolific", d_mesh, pro),
      d_nut(this, "Nutrient", d_mesh, nut), d_hyp(this, "Hypoxic", d_mesh, hyp),
      d_taf(this, "TAF", d_mesh, taf), d_ecm(this, "ECM", d_mesh, ecm),
      d_mde(this, "MDE", d_mesh, mde), d_pres(this, "Pressure", d_mesh, pres),
      d_grad_taf(this, "TAF_Gradient", d_mesh, grad_taf),
      d_vel(this, "Velocity", d_mesh, vel), d_tum(this, "Tumor", d_mesh, tum),
      d_ghosting_fv(d_mesh) {

  d_sys_names.clear();
  for (auto s : get_all_assembly()) {
    if (d_sys_names.size() <= s->d_sys.number())
      d_sys_names.resize(s->d_sys.number() + 3);
    d_sys_names[s->d_sys.number() + 2] = s->d_sys_name;
  }
  d_sys_names[0] = "Pressure_1D";
  d_sys_names[1] = "Nutrient_1D";

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
    for (auto &s : get_all_assembly())
      s->d_sys.attach_assemble_object(*s);

    // Initialize and print system
    if (!d_input.d_restart)
      d_tum_sys.init();

    for (auto &s : get_all_assembly())
      s->d_sys.time = d_input.d_init_time;

    if (d_input.d_perform_output and !d_input.d_quiet) {
      d_delayed_msg += "Libmesh Info \n";
      d_delayed_msg += d_mesh.get_info();
      d_delayed_msg += " \n";
      d_delayed_msg += d_tum_sys.get_info();
      d_mesh.write("mesh_" + d_input.d_outfile_tag + ".e");
    }

    // allocate memory to compute error for convergence
    d_err_check_pro = d_pro.d_sys.solution->clone();
    d_err_check_hyp = d_hyp.d_sys.solution->clone();
    d_err_check_nec = d_nec.d_sys.solution->clone();
    d_err_check_nut = d_nut.d_sys.solution->clone();
    d_err_check_ecm = d_ecm.d_sys.solution->clone();
    d_err_check_mde = d_mde.d_sys.solution->clone();
  }

  // 1-D network
  oss << "[Network] \n";
  d_log(oss, "init");
  d_network.create_initial_network();
  d_log(d_delayed_msg, "debug");

  // we require pressure, nutrient, and TAF localized to each processor
  d_pres.init_localized_sol(*d_comm_p);
  d_nut.init_localized_sol(*d_comm_p);
  d_taf.init_localized_sol(*d_comm_p);

  // save setup end time
  clock_end = steady_clock::now();
  d_log.d_setup_time = util::TimePair(clock_begin, clock_end);
}

void netfcfvfe::Model::run() {

  // print initial state
  if (d_input.d_perform_output) {
    d_is_output_step = true;
    write_system(0);
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

  if (!d_input.d_test_name.empty())
    libmesh_warning(
      "Test of sub-system removed from the model.\n Name of "
      "test will be discarded and instead full system will be solved.");

  do {

    // Prepare time step
    d_step++;
    d_time += d_dt;
    d_network.d_update_number += 1;

    // check if this is output step
    d_is_output_step = false;
    if ((d_step >= d_input.d_dt_output_interval and
         (d_step % d_input.d_dt_output_interval == 0)) or
        d_step == 0)
      d_is_output_step = true;

    d_is_growth_step = false;
    if (d_step % d_input.d_network_update_interval == 0 and d_input.d_network_update)
      d_is_growth_step = true;

    // init ts log
    d_log.ready_new_step(int(d_step) - 1, d_time);
    solve_clock = steady_clock::now();

    d_log("Time step: " + std::to_string(d_step) +
            ", time: " + std::to_string(d_time) + "\n",
          "integrate");

    // solve tumor-network system
    solve_system();

    // update network
    d_network.updateNetwork(d_taf, d_grad_taf);

    // write tumor solution
    write_system((d_step - d_input.d_init_step) /
                 d_input.d_dt_output_interval);
  } while (d_step < d_input.d_steps);
}

void netfcfvfe::Model::write_system(const unsigned int &t_step) {

  if (d_step >= 1) {
    // compute qoi
    compute_qoi();
    d_log.log_qoi(d_time, d_qoi.get_last());

    // add to log
    d_log.add_solve_time(util::TimePair(solve_clock, steady_clock::now()));

    // output time logger info
    d_log.log_ts();
  }

  // check if for writing to vtk files
  if (!d_is_output_step)
    return;

  // write network simulation
  d_network.writeDataToVTKTimeStep_VGM(t_step);
  d_networkVtkWriter.write( d_network.VGM, t_step );

  // write tumor simulation
  rw::VTKIO(d_mesh).write_equation_systems(
    d_input.d_outfilename + "_" + std::to_string(t_step) + ".pvtu",
    d_tum_sys);
}

void netfcfvfe::Model::compute_qoi() {

  // initialize qoi data
  int N = 17;
  std::vector<double> qoi(N, 0.);
  if (d_qoi.d_vec.empty()) {
    d_qoi = util::QoIVec(
      {"tumor_mass", "hypoxic_mass", "necrotic_mass", "prolific_mass",
       "nutrient_mass", "tumor_l2", "hypoxic_l2", "necrotic_l2",
       "prolific_l2", "nutrient_l2", "r_v_mean", "r_v_std", "l_v_mean",
       "l_v_std", "l_v_total", "vessel_vol", "vessel_density"});

    d_log.log_qoi_header(d_time, d_qoi.get_names());
  }

  qoi[0] = d_tum.compute_qoi("mass");
  qoi[1] = d_hyp.compute_qoi("mass");
  qoi[2] = d_nec.compute_qoi("mass");
  qoi[3] = d_pro.compute_qoi("mass");
  qoi[4] = d_nut.compute_qoi("mass");

  qoi[5] = d_tum.compute_qoi("l2");
  qoi[6] = d_hyp.compute_qoi("l2");
  qoi[7] = d_nec.compute_qoi("l2");
  qoi[8] = d_pro.compute_qoi("l2");
  qoi[9] = d_nut.compute_qoi("l2");

  // get other qoi such as radius mean, radius std dev,
  // length mean, length std dev, total length, total vessel vol, total domain vol
  auto vessel_qoi = d_network.compute_qoi();
  unsigned int i = 10; // start of network qoi
  for (const auto &q : vessel_qoi)
    qoi[i++] = q;

  d_qoi.add(qoi);

  // write to files for Tobias analysis
  {
    auto file_tags = d_qoi.get_names();

    for (int i = 0; i < file_tags.size(); i++) {
      std::fstream file_qoi;
      file_qoi.open("two_vessels_total_" + file_tags[i] + ".txt", std::ios::out | ios::app);
      file_qoi << qoi[i] << std::endl;
    }
  }
}

void netfcfvfe::Model::solve_system() {

  // update time
  for (auto &s : get_all_assembly())
    s->d_sys.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  for (auto &s : get_all_assembly()) {
    if (s->d_sys_name != "Velocity" and s->d_sys_name != "Tumor" and s->d_sys_name != "TAF_Gradient") {

      *(s->d_sys.old_local_solution) = *(s->d_sys.current_local_solution);
    }
  }

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // reset nonlinear step
  d_nonlinear_step = 0;

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

    // solver for 1D + 3D nutrient
    reset_clock();
    d_log("|Nutrient_1D + " + d_nut.d_sys_name + "| -> ", "solve sys");
    d_network.solve3D1DNutrientProblem(d_nut, d_tum);
    d_log.add_sys_solve_time(clock_begin, 1);
    d_log.add_sys_solve_time(clock_begin, d_nut.d_sys.number());

    reset_clock();
    d_log("|" + d_pro.d_sys_name + "| -> ", "solve sys");
    d_err_check_pro->zero();
    d_err_check_pro->add(*(d_pro.d_sys.solution));
    d_pro.solve();
    d_err_check_pro->add(-1., *(d_pro.d_sys.solution));
    d_log.add_sys_solve_time(clock_begin, d_pro.d_sys.number());

    reset_clock();
    d_log("|" + d_hyp.d_sys_name + "| -> ", "solve sys");
    d_err_check_hyp->zero();
    d_err_check_hyp->add(*(d_hyp.d_sys.solution));
    d_hyp.solve();
    d_err_check_hyp->add(-1., *(d_hyp.d_sys.solution));
    d_log.add_sys_solve_time(clock_begin, d_hyp.d_sys.number());

    reset_clock();
    d_log("|" + d_nec.d_sys_name + "| -> ", "solve sys");
    d_err_check_nec->zero();
    d_err_check_nec->add(*(d_nec.d_sys.solution));
    d_nec.solve();
    d_err_check_nec->add(-1., *(d_nec.d_sys.solution));
    d_log.add_sys_solve_time(clock_begin, d_nec.d_sys.number());

    reset_clock();
    d_log("|" + d_mde.d_sys_name + "| -> ", "solve sys");
    d_err_check_mde->zero();
    d_err_check_mde->add(*(d_mde.d_sys.solution));
    d_mde.solve();
    d_err_check_mde->add(-1., *(d_mde.d_sys.solution));
    d_log.add_sys_solve_time(clock_begin, d_mde.d_sys.number());

    reset_clock();
    d_log("|" + d_ecm.d_sys_name + "| \n", "solve sys");
    d_err_check_ecm->zero();
    d_err_check_ecm->add(*(d_ecm.d_sys.solution));
    d_ecm.solve();
    d_err_check_ecm->add(-1., *(d_ecm.d_sys.solution));
    d_log.add_sys_solve_time(clock_begin, d_ecm.d_sys.number());

    // Nonlinear iteration error
    double nonlinear_iter_error = d_err_check_pro->linfty_norm() + d_err_check_hyp->linfty_norm() + d_err_check_nec->linfty_norm() + d_err_check_mde->linfty_norm() + d_err_check_ecm->linfty_norm();

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = d_pro.d_sys.n_linear_iterations();
      const Real final_linear_residual = d_pro.d_sys.final_linear_residual();

      oss << "      LC step: " << n_linear_iterations
          << ", res: " << final_linear_residual
          << ", NC: ||u - u_old|| = "
          << nonlinear_iter_error << std::endl
          << std::endl;
      d_log(oss, "debug");
    }
    if (nonlinear_iter_error < d_input.d_nonlin_tol) {

      d_log(" \n", "debug");
      oss << "      NC step: " << l << std::endl
          << std::endl;
      d_log(oss, "converge sys");

      break;
    }
  } // nonlinear solver loop

  d_log(" \n", "solve sys");
  d_log.add_nonlin_iter(d_nonlinear_step);

  //if (d_is_growth_step or d_is_output_step) {
    // solve taf
    reset_clock();
    d_log("      Solving |" + d_taf.d_sys_name + "| \n", "solve sys");
    d_taf.solve();
    d_log.add_sys_solve_time(clock_begin, d_taf.d_sys.number());

    // Note: Grad TAF is not really used in growth algorithm
    // so we disable it
    if (false) {
      // solve for grad taf
      reset_clock();
      d_log("      Solving |" + d_grad_taf.d_sys_name + "| \n", "solve sys");
      d_grad_taf.solve();
      d_log.add_sys_solve_time(clock_begin, d_grad_taf.d_sys.number());
    }

    // solve for tumor
    d_log("      Solving |" + d_tum.d_sys_name + "| \n", "solve sys");
    d_tum.solve_custom();
    d_log.add_sys_solve_time(clock_begin, d_tum.d_sys.number());
  //}

  d_log(" \n", "solve sys");
}

void netfcfvfe::Model::solve_pressure() {

  auto solve_clock = steady_clock::now();
  reset_clock();

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_pres(
    d_pres.d_sys.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
    d_input.d_linear_tol;

  // call fully coupled 1D-3D solver (this will update the 3D pressure system
  // in libmesh
  d_log("      Solving |Pressure_1D + " + d_pres.d_sys_name + "| \n", "solve pres");
  d_network.solve3D1DFlowProblem(d_pres, d_tum);
  clock_end = std::chrono::steady_clock::now();
  if (d_log.d_cur_step >= 0) {
    d_log.add_pres_solve_time(util::TimePair(solve_clock, clock_end));
    d_log.add_pres_nonlin_iter(1);
  }

  // Note: Solving for velocity takes unusually long time so we disable it
  //  if the advection effects are disabled
  if (d_input.d_advection_active) {
    // solve for velocity
    reset_clock();
    d_log("      Solving |" + d_vel.d_sys_name + "| \n", "solve pres");
    d_log(" \n", "solve pres");
    d_vel.solve();
    if (d_log.d_cur_step >= 0)
      d_log.add_sys_solve_time(clock_begin, d_vel.d_sys.number());
  }
}
