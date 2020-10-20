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

namespace netpresnut {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(netpresnut::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Pressure")
    sys.project_solution(netpresnut::initial_condition_pres, nullptr,
                         es.parameters);
  else {
    return;
  }
}
} // namespace netpresnut

// Model setup and run
void netpresnut::model_setup_run(int argc, char **argv,
                                 const std::string &filename,
                                 Parallel::Communicator *comm) {

  auto sim_begin = steady_clock::now();

  reset_clock();

  // init seed for random number
  //random_init();

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


  log("Model: NetPresNut\n", "init");

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
  auto &pres = tum_sys.add_system<TransientLinearImplicitSystem>("Pressure");

  // some initial setups
  {
    if (input.d_restart) {
      nut.update();
      pres.update();
    }

    if (!input.d_restart) {
      // variable in nutrient system
      nut.add_variable("nutrient", CONSTANT, MONOMIAL);

      // variable in Pressure system
      pres.add_variable("pressure", CONSTANT, MONOMIAL);

      // attach initial condition function to systems
      nut.attach_init_function(initial_condition);
      pres.attach_init_function(initial_condition);
    }

    // Add boundary condition
    boundary_condition_nut(tum_sys);
  }

  //
  // Create Model class
  //
  auto model =
    Model(argc, argv, filename, comm, input, mesh, tum_sys, nut, pres, log);

  // run model
  model.run();

  model.d_log.log_ts_final(
    util::time_diff(sim_begin, steady_clock::now()));
}

// Model class
netpresnut::Model::Model(
  int argc, char **argv, const std::string &filename, Parallel::Communicator *comm,
  InpDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
  TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &pres,
  util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log, "NetPresNut"),
      d_network(this),
      d_networkVtkWriter(comm, input.d_outfilename_net + "new_"),
      d_networkVtkWriterOld(comm, input.d_outfilename_net + "old_"),
      d_nut_assembly(this, "Nutrient", d_mesh, nut),
      d_pres_assembly(this, "Pressure", d_mesh, pres),
      d_ghosting_fv(d_mesh) {

  d_nut_id = d_nut_assembly.d_sys.number();
  d_pres_id = d_pres_assembly.d_sys.number();
  d_pres_1d_id = d_pres_id + 1;
  d_nut_1d_id = d_pres_id + 2;

  d_sys_names = {"Nutrient", "Pressure", "Pressure_1D", "Nutrient_1D"};

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
    nut.attach_assemble_object(d_nut_assembly);
    pres.attach_assemble_object(d_pres_assembly);

    // add ghosting functors
    //pres.get_dof_map().add_coupling_functor(d_ghosting_fv);
    //nut.get_dof_map().add_coupling_functor(d_ghosting_fv);

    //
    // Initialize and print system
    //
    if (!d_input.d_restart)
      d_tum_sys.init();

    nut.time = d_input.d_init_time;
    pres.time = d_input.d_init_time;

    if (d_input.d_perform_output and !d_input.d_quiet) {
      d_delayed_msg += "Libmesh Info \n";
      d_delayed_msg += d_mesh.get_info();
      d_delayed_msg += " \n";
      d_delayed_msg += d_tum_sys.get_info();
      d_mesh.write("mesh_" + d_input.d_outfile_tag + ".e");
    }
  }

  // 1-D network
  oss << "[Network] \n";
  d_log(oss, "init");
  d_network.create_initial_network();
  d_log(d_delayed_msg, "debug");

  // we require pressure, nutrient, and TAF localized to each processor
  d_pres_assembly.init_localized_sol(*d_comm_p);
  d_nut_assembly.init_localized_sol(*d_comm_p);

  // save setup end time
  clock_end = steady_clock::now();
  d_log.d_setup_time = util::TimePair(clock_begin, clock_end);
}

void netpresnut::Model::run() {

  // print initial state
  {
    // Tumor model
    if (d_input.d_perform_output)
      write_system(0);

    // network
    d_networkVtkWriterOld.write(d_network.VGM, 0);
    d_networkVtkWriter.write(d_network.VGM, 0);
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

  // check for pressure test
  if (d_input.d_test_name == "test_pressure") {

    // solve for pressure only once
    auto nt = d_input.d_nonlin_max_iters;
    d_input.d_nonlin_max_iters = 2 * nt;
    solve_pressure();
    d_input.d_nonlin_max_iters = nt;

    // write tumor solution
    write_system(1);
    d_networkVtkWriterOld.write(d_network.VGM, 1.);
    d_networkVtkWriter.write(d_network.VGM, 1);

    return;
  }

  do {

    // Prepare time step
    d_step++;
    d_time += d_dt;

    // init ts log
    d_log.ready_new_step(int(d_step) - 1, d_time);
    auto solve_clock = steady_clock::now();

    // check if this is output step
    d_is_output_step = false;
    if ((d_step >= d_input.d_dt_output_interval and
         (d_step % d_input.d_dt_output_interval == 0)) or
        d_step == 0)
      d_is_output_step = true;

    d_network.d_update_number += 1;

    oss << "Time step: " << d_step << ", time: " << d_time << "\n";
    d_log(oss, "integrate");
    d_log(" \n", "integrate");

    // solve tumor-network system
    solve_system();

    // compute qoi
    compute_qoi();

    // Post-processing
    if (d_is_output_step) {

      // write tumor solution
      write_system((d_step - d_input.d_init_step) /
                   d_input.d_dt_output_interval);

      const int timeStep = static_cast<int>((d_step - d_input.d_init_step) / d_input.d_dt_output_interval);

      d_networkVtkWriterOld.write(d_network.VGM, timeStep);
      d_networkVtkWriter.write(d_network.VGM, timeStep);
    }

    // output qoi
    if (d_step == 1)
      d_log.log_qoi_header(d_time, d_qoi.get_names());
    d_log.log_qoi(d_time, d_qoi.get_last());

    // add to log
    d_log.add_solve_time(util::TimePair(solve_clock, steady_clock::now()));

    // output time logger info
    d_log.log_ts();

  } while (d_step < d_input.d_steps);
}

void netpresnut::Model::write_system(const unsigned int &t_step) {

  //
  rw::VTKIO(d_mesh).write_equation_systems(
    d_input.d_outfilename + "_" + std::to_string(t_step) + ".pvtu",
    d_tum_sys);
}

void netpresnut::Model::compute_qoi() {

  return;
}

void netpresnut::Model::solve_system() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &pres = d_tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");

  // update time
  nut.time = d_time;
  pres.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *pres.old_local_solution = *pres.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // reset nonlinear step
  d_nonlinear_step = 0;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_nut(
    nut.solution->clone());

  d_log("  Nonlinear loop\n", "solve sys");

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
    d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    Nonlinear step: " << l << " ";
    d_log(oss, "solve sys");

    // solver for 1D nutrient
    d_log("|1D nutrient| -> ", "solve sys");
    d_network.solveVGMforNutrient(d_pres_assembly, d_nut_assembly);
    d_log.add_sys_solve_time(clock_begin, d_nut_1d_id);

    // solve nutrient
    reset_clock();
    last_nonlinear_soln_nut->zero();
    last_nonlinear_soln_nut->add(*nut.solution);
    d_log("|3D nutrient| \n", "solve sys");
    nut.solve();
    last_nonlinear_soln_nut->add(-1., *nut.solution);
    last_nonlinear_soln_nut->close();
    d_log.add_sys_solve_time(clock_begin, d_nut_id);

    // Nonlinear iteration error
    double nonlinear_iter_error = last_nonlinear_soln_nut->linfty_norm();
    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = nut.n_linear_iterations();
      const Real final_linear_residual = nut.final_linear_residual();

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
}

void netpresnut::Model::solve_pressure() {

  auto solve_clock = steady_clock::now();
  reset_clock();

  // get systems
  auto &pres = d_tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");

  // update time
  pres.time = d_time;

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
    d_log(" |1D pressure| -> ", "solve pres");
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
    double nonlinear_iter_error = last_nonlinear_soln_pres->linfty_norm();
    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = pres.n_linear_iterations();
      const Real final_linear_residual = pres.final_linear_residual();

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

  d_log(" \n", "solve pres");

  clock_end = std::chrono::steady_clock::now();
  if (d_log.d_cur_step >= 0) {
    d_log.add_pres_solve_time(util::TimePair(solve_clock, clock_end));
    d_log.add_pres_nonlin_iter(d_nonlinear_step);
  }
}