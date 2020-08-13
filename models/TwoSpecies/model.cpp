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

namespace twosp {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(twosp::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Tumor")
    sys.project_solution(twosp::initial_condition_tum, nullptr, es.parameters);
  else {
    return;
  }
}
} // namespace twosp

// Model setup and run
void twosp::model_setup_run(int argc, char **argv,
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


  log("Model: TwoSpecies\n", "init");

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

  // some initial setups
  {
    if (input.d_restart) {
      nut.update();
      tum.update();
    }

    if (!input.d_restart) {
      // variable in nutrient system
      nut.add_variable("nutrient", FIRST);

      // variable in tumor system
      tum.add_variable("tumor", FIRST);
      tum.add_variable("chemical_tumor", FIRST);

      // attach initial condition function to systems
      tum.attach_init_function(initial_condition);
      nut.attach_init_function(initial_condition);
    }

    // Add boundary condition
    boundary_condition_nut(tum_sys);
  }

  //
  // Create Model class
  //
  auto model = Model(argc, argv, filename, comm, input, mesh, tum_sys,
                     tum, nut, log);

  // run model
  model.run();

  model.d_log.log_ts_final(
      util::time_diff(sim_begin, steady_clock::now()));
}

// Model class
twosp::Model::Model(
    int argc, char **argv, const std::string &filename, Parallel::Communicator *comm,
    InpDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
    TransientLinearImplicitSystem &tum, TransientLinearImplicitSystem &nut,
    util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log, "TwoSpecies"),
      d_tum_assembly(this, "Tumor", d_mesh, tum),
      d_nut_assembly(this, "Nutrient", d_mesh, nut) {

  d_nut_id = d_nut_assembly.d_sys.number();
  d_tum_id = d_tum_assembly.d_sys.number();

  d_sys_names.resize(2);
  d_sys_names[d_nut_id] = "Nutrient";
  d_sys_names[d_tum_id] = "Tummor";

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
    tum.attach_assemble_object(d_tum_assembly);
    nut.attach_assemble_object(d_nut_assembly);

    //
    // Initialize and print system
    //
    if (!d_input.d_restart)
      d_tum_sys.init();

    nut.time = d_input.d_init_time;
    tum.time = d_input.d_init_time;

    if (d_input.d_perform_output and !d_input.d_quiet) {
      d_delayed_msg += "Libmesh Info \n";
      d_delayed_msg += d_mesh.get_info();
      d_delayed_msg += " \n";
      d_delayed_msg +=d_tum_sys.get_info();
      d_mesh.write("mesh_" + d_input.d_outfile_tag + ".e");
    }
  }

  // save setup end time
  clock_end = steady_clock::now();
  d_log.d_setup_time = util::TimePair(clock_begin, clock_end);
}

void twosp::Model::run() {

  // print initial state
  {
    // Tumor model
    if (d_input.d_perform_output)
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

void twosp::Model::write_system(const unsigned int &t_step) {

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
}

void twosp::Model::compute_qoi() {

  // initialize qoi data
  int N = 3;
  std::vector<double> qoi(N, 0.);
  if (d_step <= 1)
    d_qoi = util::QoIVec({"tumor_mass", "nutrient_mass", "tumor_l2"});

  // integral of total tumor
  double value_mass = 0.;
  double total_mass = 0.;
  util::computeMass(d_tum_sys, "Tumor", "tumor", value_mass);
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  qoi[0] = total_mass;

  // integral of hypoxic
  value_mass = 0.; total_mass = 0.;
  util::computeMass(d_tum_sys, "Nutrient", "nutrient", value_mass);
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  qoi[1] = total_mass;

  // L2 norm of total tumor
  value_mass = 0.; total_mass = 0.;
  value_mass = std::pow(d_input.d_mesh_size, 1.5) * d_tum_assembly.d_sys.solution->l2_norm();
  //MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE,
  //              MPI_SUM, MPI_COMM_WORLD);
  qoi[3] = value_mass;

  d_qoi.add(qoi);
}

void twosp::Model::solve_system() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");

  // update time
  nut.time = d_time;
  tum.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;

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
    d_log("|nutrient| -> ", "solve sys");
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