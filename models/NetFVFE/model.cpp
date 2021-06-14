///////////////////////////////////////////////////////////////////////////////
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

void reset_clock() { clock_begin = steady_clock::now(); }

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
    sys.project_solution(netfvfe::initial_condition_nut, nullptr,
                         es.parameters);
  else if (system_name == "Prolific")
    sys.project_solution(netfvfe::initial_condition_pro, nullptr,
                         es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(netfvfe::initial_condition_hyp, nullptr,
                         es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(netfvfe::initial_condition_nec, nullptr,
                         es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(netfvfe::initial_condition_taf, nullptr,
                         es.parameters);
  else if (system_name == "ECM")
    sys.project_solution(netfvfe::initial_condition_ecm, nullptr,
                         es.parameters);
  else if (system_name == "MDE")
    sys.project_solution(netfvfe::initial_condition_mde, nullptr,
                         es.parameters);
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

  auto sim_begin = steady_clock::now();

  reset_clock();

  // init seed for random number
  // random_init();

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
  auto &pro = tum_sys.add_system<TransientLinearImplicitSystem>("Prolific");
  auto &hyp = tum_sys.add_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = tum_sys.add_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = tum_sys.add_system<TransientLinearImplicitSystem>("TAF");
  auto &ecm = tum_sys.add_system<TransientLinearImplicitSystem>("ECM");
  auto &mde = tum_sys.add_system<TransientLinearImplicitSystem>("MDE");
  auto &pres = tum_sys.add_system<TransientLinearImplicitSystem>("Pressure");
  auto &grad_taf =
    tum_sys.add_system<TransientLinearImplicitSystem>("TAF_Gradient");
  auto &vel = tum_sys.add_system<TransientLinearImplicitSystem>("Velocity");
  auto &tum = tum_sys.add_system<TransientLinearImplicitSystem>("Tumor");

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
  auto model = Model(argc, argv, filename, comm, input, mesh, tum_sys, nec, pro,
                     nut, hyp, taf, ecm, mde, pres, grad_taf, vel, tum, log);

  // run model
  model.run();

  model.d_log.log_ts_final(util::time_diff(sim_begin, steady_clock::now()));
}

// Model class
netfvfe::Model::Model(
  int argc, char **argv, const std::string &filename,
  Parallel::Communicator *comm, InpDeck &input, ReplicatedMesh &mesh,
  EquationSystems &tum_sys, TransientLinearImplicitSystem &nec,
  TransientLinearImplicitSystem &pro, TransientLinearImplicitSystem &nut,
  TransientLinearImplicitSystem &hyp, TransientLinearImplicitSystem &taf,
  TransientLinearImplicitSystem &ecm, TransientLinearImplicitSystem &mde,
  TransientLinearImplicitSystem &pres,
  TransientLinearImplicitSystem &grad_taf, TransientLinearImplicitSystem &vel,
  TransientLinearImplicitSystem &tum, util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log, "NetFVFE"),
      d_network(this),
      d_networkVtkWriter(comm, input.d_outfilename_net + ""),
      d_networkDGFWriter(comm, input.d_outfilename_net + ""),
      d_qoi_writer(comm, "qoi.m"),
      d_csv_qoi_writer(comm, "qoi.csv"),
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
    d_err_check_pres = d_pres.d_sys.solution->clone();
  }

  // 1-D network
  oss << "[Network] \n";
  d_log(oss, "init");
  d_network.create_initial_network();
  d_network.d_extrapolate_nutrients_at_tips = d_input.d_extrapolate_nutrients_at_tips;
  d_log(d_delayed_msg, "debug");

  // we require pressure, nutrient, and TAF localized to each processor
  d_pres.init_localized_sol(*d_comm_p);
  d_nut.init_localized_sol(*d_comm_p);
  d_taf.init_localized_sol(*d_comm_p);

  // for velocity system, assembly matrix only once and then re-use
  // this restricts libmesh from zeroing matrix and vector
  d_vel.d_sys.zero_out_matrix_and_rhs = false;

  // save setup end time
  clock_end = steady_clock::now();
  d_log.d_setup_time = util::TimePair(clock_begin, clock_end);
}

void netfvfe::Model::run() {

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

  // find solver type and set assembly type for each system
  if (d_input.d_scheme_name == "solve_explicit") {
    d_nut.d_implicit_assembly = false;
    d_pro.d_implicit_assembly = false;
    d_hyp.d_implicit_assembly = false;
    d_nec.d_implicit_assembly = false;
    d_ecm.d_implicit_assembly = false;
    d_mde.d_implicit_assembly = false;
    d_taf.d_implicit_assembly = false;
  } else if (d_input.d_scheme_name == "solve_nutrient_explicit") {
    d_nut.d_implicit_assembly = false;
  }

  //
  // Solve step
  //

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
    if (d_step % d_input.d_network_update_interval == 0 and
        d_input.d_network_update)
      d_is_growth_step = true;

    // init ts log
    d_log.ready_new_step(int(d_step) - 1, d_time);
    solve_clock = steady_clock::now();

    d_log("Time step: " + std::to_string(d_step) +
            ", time: " + std::to_string(d_time) + "\n",
          "integrate");

    // stochastic contributions from the cylindrical Wiener process
    d_hyp.calculate_new_stochastic_coefficients(d_dt);
    d_pro.calculate_new_stochastic_coefficients(d_dt);

    // project prolific and hypoxic to physical range
    if (util::locate_in_set<std::string>("hypoxic", d_input.d_project_fields))
      d_hyp.project_physical_range();
    if (util::locate_in_set<std::string>("prolific", d_input.d_project_fields))
      d_pro.project_physical_range();
    if (util::locate_in_set<std::string>("necrotic", d_input.d_project_fields))
      d_nec.project_physical_range();
    if (util::locate_in_set<std::string>("tumor", d_input.d_project_fields))
      d_tum.project_physical_range();

    // solve tumor-network system
    solve_system();


    if (d_is_growth_step) {
      d_network.updateNetwork(d_taf, d_grad_taf);
    }

    // write tumor solution
    write_system(d_step); //(d_step - d_input.d_init_step) / d_input.d_dt_output_interval);
  } while (d_step < d_input.d_steps);
}

void netfvfe::Model::write_system(const unsigned int &t_step) {

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
  //if (!(t_step % 3 == 0)) // hard-code --> specify interval in input file using 'output_interval'
  d_is_output_step = false;
  if ((d_step >= d_input.d_dt_output_interval and
       (d_step % d_input.d_dt_output_interval == 0)) or
      d_step == 0)
    d_is_output_step = true;

  if (!d_is_output_step)
    return;

  unsigned int file_number = d_step; // using actual time step
  // int file_number = (d_step - d_input.d_init_step) / d_input.d_dt_output_interval;

  // write network simulation
  d_networkVtkWriter.write(d_network.VGM, file_number);
  d_networkDGFWriter.write(d_network.VGM);

  // write tumor simulation
  rw::VTKIO(d_mesh).write_equation_systems(d_input.d_outfilename + "_" +
                                             std::to_string(file_number) + ".pvtu",
                                           d_tum_sys);
}

void netfvfe::Model::compute_qoi() {

  std::vector<std::string> qoi_names =  {"tumor_mass", "hypoxic_mass", "necrotic_mass", "prolific_mass",
                                         "nutrient_mass", "tumor_l2", "hypoxic_l2", "necrotic_l2",
                                         "prolific_l2", "nutrient_l2", "r_v_mean", "r_v_std", "l_v_mean",
                                         "l_v_std", "l_v_total", "vessel_vol", "vessel_density",
                                         "network_total_added_length", "network_total_removed_length",
                                         "network_total_added_volume", "network_total_removed_volume",
                                         "network_total_length", "network_total_volume",
                                         "min_tumor", "max_tumor",
                                         "min_hypoxic", "max_hypoxic",
                                         "min_prolific", "max_prolific",
                                         "min_necrotic", "max_necrotic",
                                         "time" };

  // initialize qoi data
  const auto N = qoi_names.size();
  std::vector<double> qoi(N, 0.);
  if (d_qoi.d_vec.empty()) {
    d_qoi = util::QoIVec(qoi_names);

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
  // length mean, length std dev, total length, total vessel vol, total domain
  // vol
  auto vessel_qoi = d_network.compute_qoi();
  unsigned int i = 10; // start of network qoi
  for (const auto &q : vessel_qoi)
    qoi[i++] = q;

  // if it was a growth step we add statistical data how much the network changed
  if (d_is_growth_step) {
    qoi[i + 0] = d_network.total_added_length;
    qoi[i + 1] = d_network.total_removed_length;
    qoi[i + 2] = d_network.total_added_volume;
    qoi[i + 3] = d_network.total_removed_volume;
  } else {
    qoi[i + 0] = 0;
    qoi[i + 1] = 0;
    qoi[i + 2] = 0;
    qoi[i + 3] = 0;
  }

  // we add the length and volume of the 1D network
  double length{0}, volume{0};
  d_network.get_length_and_volume_of_network(length, volume);
  qoi[i + 4] = length;
  qoi[i + 5] = volume;

  // we save the min and max values of the tumor, hypoxic and prolific concentrations.
  i = i + 6;
  {
    auto min = d_tum.d_sys.solution->min();
    d_mesh.comm().min(min);
    auto max = d_tum.d_sys.solution->max();
    d_mesh.comm().max(max);
    qoi[i++] = min;
    qoi[i++] = max;
  }

  {
    auto min = d_hyp.d_sys.solution->min();
    d_mesh.comm().min(min);
    auto max = d_hyp.d_sys.solution->max();
    d_mesh.comm().max(max);
    qoi[i++] = min;
    qoi[i++] = max;
  }

  {
    auto min = d_pro.d_sys.solution->min();
    d_mesh.comm().min(min);
    auto max = d_pro.d_sys.solution->max();
    d_mesh.comm().max(max);
    qoi[i++] = min;
    qoi[i++] = max;
  }

  {
    auto min = d_nec.d_sys.solution->min();
    d_mesh.comm().min(min);
    auto max = d_nec.d_sys.solution->max();
    d_mesh.comm().max(max);
    qoi[i++] = min;
    qoi[i++] = max;
  }

  qoi[i++] = d_time;

  d_qoi.add(qoi);

  d_qoi_writer.write(d_qoi);
  d_csv_qoi_writer.write(d_qoi);
}

void netfvfe::Model::solve_system() {

  if (d_input.d_scheme_name == "solve_explicit")
    solve_system_explicit();
  else if (d_input.d_scheme_name == "solve_nutrient_explicit")
    solve_system_nutrient_explicit();
  else if (d_input.d_scheme_name == "solve_implicit")
    solve_system_implicit();
  else
    libmesh_error_msg("Scheme name is not valid. Try solve_explicit, "
                      "solve_nutrient_explicit or solve_implicit.");
}

void netfvfe::Model::solve_system_explicit() {

  // update time
  for (auto &s : get_all_assembly())
    s->d_sys.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  std::vector<std::string> no_update_sys = {"Pressure", "Velocity", "Tumor",
                                            "TAF_Gradient"};
  if (!d_input.d_solve_ecm) {
    no_update_sys.push_back("MDE");
    no_update_sys.push_back("ECM");
  }

  for (auto &s : get_all_assembly()) {
    if (util::locate_in_set(s->d_sys_name, no_update_sys) == -1)
      *(s->d_sys.old_local_solution) = *(s->d_sys.current_local_solution);
  }

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // solve nutrient
  solve_nutrient();

  reset_clock();
  d_log("Solving |" + d_pro.d_sys_name + "| -> ", "solve sys");
  d_pro.solve();
  d_log.add_sys_solve_time(clock_begin, d_pro.d_sys.number());

  reset_clock();
  d_log("|" + d_hyp.d_sys_name + "| -> ", "solve sys");
  d_hyp.solve();
  d_log.add_sys_solve_time(clock_begin, d_hyp.d_sys.number());

  reset_clock();
  d_log("|" + d_nec.d_sys_name + "| -> ", "solve sys");
  d_nec.solve();
  d_log.add_sys_solve_time(clock_begin, d_nec.d_sys.number());

  // solve mde and ecm only if specified in the input file
  if (d_input.d_solve_ecm) {
    reset_clock();
    d_log("|" + d_mde.d_sys_name + "| -> ", "solve sys");
    d_mde.solve();
    d_log.add_sys_solve_time(clock_begin, d_mde.d_sys.number());

    reset_clock();
    d_log("|" + d_ecm.d_sys_name + "| -> ", "solve sys");
    d_ecm.solve();
    d_log.add_sys_solve_time(clock_begin, d_ecm.d_sys.number());
  }

  // solve taf and tumor only if it is either growth step or output step
  //if (d_is_growth_step or d_is_output_step) {
  // solve taf
  reset_clock();
  d_log("|" + d_taf.d_sys_name + "| -> ", "solve sys");
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
  d_log("|" + d_tum.d_sys_name + "| \n", "solve sys");
  d_tum.solve_custom();
  d_log.add_sys_solve_time(clock_begin, d_tum.d_sys.number());
  //}

  d_log(" \n", "solve sys");
}

void netfvfe::Model::solve_system_implicit() {

  // update time
  for (auto &s : get_all_assembly())
    s->d_sys.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  std::vector<std::string> no_update_sys = {"Pressure", "Velocity", "Tumor",
                                            "TAF_Gradient"};
  if (!d_input.d_solve_ecm) {
    no_update_sys.push_back("MDE");
    no_update_sys.push_back("ECM");
  }

  for (auto &s : get_all_assembly()) {
    if (util::locate_in_set(s->d_sys_name, no_update_sys) == -1)
      *(s->d_sys.old_local_solution) = *(s->d_sys.current_local_solution);
  }

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // tissue systems to solve
  std::vector<util::BaseAssembly *> sys_solve = {&d_pro, &d_hyp, &d_nec};
  if (d_input.d_solve_ecm) {
    sys_solve.push_back(&d_mde);
    sys_solve.push_back(&d_ecm);
  }

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

    if (d_input.d_coupled_1d3d) {

      // solver for 1D + 3D nutrient
      reset_clock();
      d_log(" |Nutrient_1D + " + d_nut.d_sys_name + "| -> ", "solve sys");
      d_network.solve3D1DNutrientProblem(d_nut, d_tum);
      d_log.add_sys_solve_time(clock_begin, 1);
      d_log.add_sys_solve_time(clock_begin, d_nut.d_sys.number());
    } else {

      // solver for nutrients
      reset_clock();
      d_log("|Nutrient_1D| -> ", "solve sys");
      d_network.solveVGMforNutrient(d_pres, d_nut);
      d_log.add_sys_solve_time(clock_begin, 1);

      reset_clock();
      d_log(" |" + d_nut.d_sys_name + "| -> ", "solve sys");
      d_err_check_nut->zero();
      d_err_check_nut->add(*(d_nut.d_sys.solution));
      d_nut.solve();
      d_err_check_nut->add(-1., *(d_nut.d_sys.solution));
      d_log.add_sys_solve_time(clock_begin, d_nut.d_sys.number());
    }

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

    if (d_input.d_solve_ecm) {
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
    }

    // Nonlinear iteration error
    double nonlinear_iter_error = d_err_check_pro->linfty_norm() + d_err_check_hyp->linfty_norm() + d_err_check_nec->linfty_norm();
    if (d_input.d_solve_ecm)
      nonlinear_iter_error += d_err_check_mde->linfty_norm() + d_err_check_ecm->linfty_norm();
    if (!d_input.d_coupled_1d3d)
      nonlinear_iter_error += d_err_check_nut->linfty_norm();

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

void netfvfe::Model::solve_system_nutrient_explicit() {

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

  // solve nutrient
  solve_nutrient();

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

    reset_clock();
    d_log(" |" + d_pro.d_sys_name + "| -> ", "solve sys");
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

    if (d_input.d_solve_ecm) {
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
    }

    // Nonlinear iteration error
    double nonlinear_iter_error = d_err_check_pro->linfty_norm() + d_err_check_hyp->linfty_norm() + d_err_check_nec->linfty_norm();
    if (d_input.d_solve_ecm)
      nonlinear_iter_error += d_err_check_mde->linfty_norm() + d_err_check_ecm->linfty_norm();

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

void netfvfe::Model::solve_pressure() {

  // disable matrix assembly for velocity
  if (d_step > 1)
    d_vel.d_assemble_matrix = false;

  // If d_solve_pres_with_net_update = true, solve for pressure only after
  // the network update
  /*
  bool solve_pres = false;
  if (d_input.d_solve_pres_with_net_update) {
    // two cases: network update is allowed or network is static
    if (d_input.d_network_update) {
      if ((d_step - 1) % d_input.d_network_update_interval == 0)
        solve_pres = true;
    } else {
      // network is static. In this case, we update pressure every 4 time step
      // FIXME: Should we introduce a input parameter for this as well??
      if (d_step % 4 == 0)
        solve_pres = true;
    }
  } else
    solve_pres = true;

  // if this is 1st time step, we have to solve for pressure in any case
  if (d_step == 1)
    solve_pres = true;

  if (!solve_pres)
    return;

  auto solve_clock = steady_clock::now();
  reset_clock();
  */

  // coupled 1d-3d or iterative method
  if (d_input.d_coupled_1d3d) {
    // call fully coupled 1D-3D solver (this will update the 3D pressure system
    // in libmesh
    d_log("      Solving |Pressure_1D + " + d_pres.d_sys_name + "| \n",
          "solve pres");
    d_network.solve3D1DFlowProblem(d_pres, d_tum);
    clock_end = std::chrono::steady_clock::now();
    if (d_log.d_cur_step >= 0) {
      d_log.add_pres_solve_time(util::TimePair(solve_clock, clock_end));
      d_log.add_pres_nonlin_iter(1);
    }
  } else {

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
      d_log(" |Pressure_1D| -> ", "solve pres");
      d_network.solveVGMforPressure(d_pres);
      if (d_log.d_cur_step >= 0)
        d_log.add_sys_solve_time(clock_begin, 0);

      // solve for pressure in tissue
      reset_clock();
      d_err_check_pres->zero();
      d_err_check_pres->add(*d_pres.d_sys.solution);
      d_log("|" + d_pres.d_sys_name + "|\n", "solve pres");
      d_pres.d_sys.solve();

      d_err_check_pres->add(-1., *d_pres.d_sys.solution);
      d_err_check_pres->close();
      if (d_log.d_cur_step >= 0)
        d_log.add_sys_solve_time(clock_begin, d_pres.d_sys.number());

      // Nonlinear iteration error
      double nonlinear_iter_error = d_err_check_pres->linfty_norm();
      if (d_input.d_perform_output) {

        const unsigned int n_linear_iterations =
          d_pres.d_sys.n_linear_iterations();
        const Real final_linear_residual = d_pres.d_sys.final_linear_residual();

        oss << "      LC step: " << n_linear_iterations
            << ", res: " << final_linear_residual
            << ", NC: ||u - u_old|| = " << nonlinear_iter_error << std::endl
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

  // Note: Solving for velocity takes unusually long time so we disable it
  //  if the advection effects are disabled
  if (d_input.d_advection_active) {
    // solve for velocity
    reset_clock();
    d_log("      Solving |" + d_vel.d_sys_name + "| \n", "solve pres");
    d_vel.solve();
    if (d_log.d_cur_step >= 0)
      d_log.add_sys_solve_time(clock_begin, d_vel.d_sys.number());
  }
}

void netfvfe::Model::solve_nutrient() {

  auto solve_clock = steady_clock::now();
  reset_clock();

  // coupled 1d-3d or iterative method
  if (d_input.d_coupled_1d3d) {
    // solver for 1D + 3D nutrient
    reset_clock();
    d_log("      Solving |Nutrient_1D + " + d_nut.d_sys_name + "| -> ", "solve nut");
    d_network.solve3D1DNutrientProblem(d_nut, d_tum);
    if (d_log.d_cur_step >= 0) {
      d_log.add_sys_solve_time(clock_begin, 1);
      d_log.add_sys_solve_time(clock_begin, d_nut.d_sys.number());
    }

  } else {
    // Nonlinear iteration loop
    d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

    // reset nonlinear step
    d_nonlinear_step = 0;

    // Solve pressure system
    d_log("  Nonlinear loop for nutrient\n\n", "solve nut");
    d_log(" \n", "solve nut");
    // nonlinear loop
    for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
      // Debug
      // for (unsigned int l = 0; l < 10; ++l) {

      d_nonlinear_step = l;
      oss << "    Nonlinear step: " << l;
      d_log(oss, "solve pres");

      // solver for nutrients
      reset_clock();
      d_log(" |Nutrient_1D| -> ", "solve nut");
      d_network.solveVGMforNutrient(d_pres, d_nut);
      if (d_log.d_cur_step >= 0)
        d_log.add_sys_solve_time(clock_begin, 1);

      reset_clock();
      d_log("|" + d_nut.d_sys_name + "| -> ", "solve nut");
      d_err_check_nut->zero();
      d_err_check_nut->add(*d_nut.d_sys.solution);
      d_nut.solve();
      d_err_check_nut->add(-1., *d_nut.d_sys.solution);
      d_err_check_nut->close();
      if (d_log.d_cur_step >= 0)
        d_log.add_sys_solve_time(clock_begin, d_nut.d_sys.number());

      // Nonlinear iteration error
      double nonlinear_iter_error = d_err_check_nut->linfty_norm();
      if (d_input.d_perform_output) {

        const unsigned int n_linear_iterations =
          d_nut.d_sys.n_linear_iterations();
        const Real final_linear_residual = d_nut.d_sys.final_linear_residual();

        oss << "      LC step: " << n_linear_iterations
            << ", res: " << final_linear_residual
            << ", NC: ||u - u_old|| = " << nonlinear_iter_error << std::endl
            << std::endl;
        d_log(oss, "debug");
      }
      if (nonlinear_iter_error < d_input.d_nonlin_tol) {

        d_log(" \n", "debug");
        oss << "      NC step: " << l << std::endl
            << std::endl;
        d_log(oss, "converge nut");

        break;
      }
    } // nonlinear solver loop

    d_log(" \n", "solve nut");
  }
}
