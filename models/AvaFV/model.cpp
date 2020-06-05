////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "model.hpp"
#include "rw/vtk_io.hpp"

namespace {
std::ostringstream oss;

void random_init() { srand(time(nullptr)); }
} // namespace

namespace avafv {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(avafv::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Tumor")
    sys.project_solution(avafv::initial_condition_tum, nullptr, es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(avafv::initial_condition_hyp, nullptr, es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(avafv::initial_condition_nec, nullptr, es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(avafv::initial_condition_taf, nullptr, es.parameters);
  else {
    return;
  }
}
} // namespace avafv

// Model setup and run
void avafv::model_setup_run(int argc, char **argv,
                             const std::string &filename,
                             Parallel::Communicator *comm) {

  // init seed for random number
  random_init();

  // test util inp
  auto input = InpDeck(filename);

  // create logger
  util::Logger log("sim_" + input.d_outfile_tag + ".log", comm, !input.d_quiet);

  // disable reference counter information
  if (input.d_quiet)
    ReferenceCounter::disable_print_counter_info();

  //
  oss << " ********** AvaFV **************\n";
  log(oss);

  // create mesh
  oss << "Creating tumor mesh\n";
  log(oss);
  ReplicatedMesh mesh(*comm);
  if (input.d_read_mesh_flag)
    mesh.read(input.d_mesh_filename);
  else
    util::create_mesh(input, mesh);

  oss << "Creating tumor system\n";
  log(oss);
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
  auto &grad_taf =
      tum_sys.add_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // some initial setups
  {
    if (input.d_restart) {
      nut.update();
      tum.update();
      hyp.update();
      nec.update();
      taf.update();
      grad_taf.update();
    }

    if (!input.d_restart) {
      // variable in nutrient system
      nut.add_variable("nutrient", CONSTANT, MONOMIAL);

      // variable in tumor system
      tum.add_variable("tumor", CONSTANT, MONOMIAL);
      tum.add_variable("chemical_tumor", CONSTANT, MONOMIAL);

      // variable in hypoxic system
      hyp.add_variable("hypoxic", CONSTANT, MONOMIAL);

      // variable in necrotic system
      nec.add_variable("necrotic", CONSTANT, MONOMIAL);

      // variable in TAF system
      taf.add_variable("taf", CONSTANT, MONOMIAL);

      // variable in gradient TAF system
      grad_taf.add_variable("taf_gradx", CONSTANT, MONOMIAL);
      grad_taf.add_variable("taf_grady", CONSTANT, MONOMIAL);
      if (input.d_dim > 2)
        grad_taf.add_variable("taf_gradz", CONSTANT, MONOMIAL);

      // attach initial condition function to systems
      tum.attach_init_function(initial_condition);
      hyp.attach_init_function(initial_condition);
      nut.attach_init_function(initial_condition);
    }

    // Add boundary condition
    boundary_condition_nut(tum_sys);
  }

  //
  // Create Model class
  //
  auto model = Model(argc, argv, filename, comm, input, mesh, tum_sys,
                     nec, tum, nut, hyp, taf, grad_taf, log);

  // run model
  model.run();
}

// Model class
avafv::Model::Model(
    int argc, char **argv, const std::string &filename,
    Parallel::Communicator *comm,
    InpDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
    TransientLinearImplicitSystem &nec, TransientLinearImplicitSystem &tum,
    TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &hyp,
    TransientLinearImplicitSystem &taf, TransientLinearImplicitSystem &grad_taf,
    util::Logger &log)
    : util::BaseModel(comm, input, mesh, tum_sys, log),
      d_nec_assembly(this, "Necrotic", d_mesh, nec),
      d_tum_assembly(this, "Tumor", d_mesh, tum),
      d_nut_assembly(this, "Nutrient", d_mesh, nut),
      d_hyp_assembly(this, "Hypoxic", d_mesh, hyp),
      d_taf_assembly(this, "TAF", d_mesh, taf),
      d_grad_taf_assembly(this, "TAF_Gradient", d_mesh, grad_taf) {

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
    // Note: For grad taf, we do not need to perform assembly
    // simply compute grad of taf at center of each element and update grad_taf
    // grad_taf.attach_assemble_object(d_grad_taf_assembly);

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
    grad_taf.time = d_input.d_init_time;

    if (d_input.d_perform_output and !d_input.d_quiet) {
      d_mesh.print_info();
      d_tum_sys.print_info();
      d_mesh.write("mesh_" + d_input.d_outfile_tag + ".e");
    }

    // set Petsc matrix option to suppress the error
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(nut.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(tum.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(hyp.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(nec.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(taf.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
  }
}

void avafv::Model::run() {

  // print initial state
  {
    //
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
  // Solve
  //

  if (!d_input.d_test_name.empty()) {
    oss << "\n\nSolving sub-system for test: " << d_input.d_test_name << "\n\n";
    d_log(oss);
  }

  // based on test_name, decide what systems to output
  if (!d_input.d_test_name.empty()) {

    d_grad_taf_assembly.d_sys.hide_output() = true;
    d_taf_assembly.d_sys.hide_output() = true;
  }

  do {

    // Prepare time step
    d_step++;
    d_time += d_dt;

    // check if this is output step
    d_is_output_step = false;
    if ((d_step >= d_input.d_dt_output_interval and
         (d_step % d_input.d_dt_output_interval == 0)) or
        d_step == 0)
      d_is_output_step = true;

    oss << "\n\n________________________________________________________\n";
    oss << "At time step: " << d_step << ", time: " << d_time << "\n\n";
    d_log(oss);

    // solve tumor-network system
    solve_system();

    // compute qoi
    compute_qoi();

    // Post-processing
    if (d_is_output_step)
      // write tumor solution
      write_system((d_step - d_input.d_init_step) /
                       d_input.d_dt_output_interval);


    // output qoi
    if (d_step == 1)
      d_log.log_qoi_header(d_time, d_qoi.get_last(), d_qoi.get_names());
    else
      d_log.log_qoi(d_time, d_qoi.get_last());

  } while (d_step < d_input.d_steps);
}

void avafv::Model::write_system(const unsigned int &t_step) {

  ExodusII_IO exodus(d_mesh);

  // write mesh and simulation results
  std::string filename = d_input.d_outfilename + ".e";

  // write to exodus
  // exodus.write_timestep(filename, d_tum_sys, 1, d_time);

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

void avafv::Model::compute_qoi() {

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

void avafv::Model::solve_system() {

  if (!d_input.d_test_name.empty()) {
    if (d_input.d_test_name == "test_tum")
      test_tum();
    else if (d_input.d_test_name == "test_tum_2")
      test_tum_2();
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
  *taf.old_local_solution = *taf.current_local_solution;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  oss << "\n  Nonlinear loop\n\n";
  d_log(oss);

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    ____________________\n";
    oss << "    Nonlinear step: " << l << "\n\n";
    oss << "      Solving ";
    d_log(oss);

    // solve nutrient
    oss << "[3D nutrient] -> ";
    d_log(oss);
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    oss << "[tumor species] -> ";
    d_log(oss);
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    oss << "[hypoxic species] -> ";
    d_log(oss);
    hyp.solve();

    // solve necrotic
    oss << "[necrotic species] -> ";
    d_log(oss);
    nec.solve();

    // solve taf
    oss << "[taf species] -> ";
    d_log(oss);
    taf.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      Linear converged at step: " << n_linear_iterations
          << ", residual: " << final_linear_residual
          << ", Nonlinear convergence: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss);
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      oss << "      Nonlinear converged at step: " << l << std::endl
          << std::endl;
      d_log(oss);

      break;
    }
  } // nonlinear solver loop

  oss << "\n  End of nonlinear loop\n";
  d_log(oss);

  // solve for gradient of taf
  oss << "      Solving [gradient of taf]\n";
  d_log(oss);
  d_grad_taf_assembly.solve();
}

void avafv::Model::test_tum() {

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

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  oss << "\n  Nonlinear loop\n";
  d_log(oss);

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    ____________________\n";
    oss << "    Nonlinear step: " << l << "\n\n";
    oss << "      Solving ";
    d_log(oss);

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    oss << "[tumor species] -> ";
    d_log(oss);
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    oss << "[hypoxic species] -> ";
    d_log(oss);
    hyp.solve();

    // solve necrotic
    oss << "[necrotic species]\n";
    d_log(oss);
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      Linear converged at step: " << n_linear_iterations
          << ", residual: " << final_linear_residual
          << ", Nonlinear convergence: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss);
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      oss << "      Nonlinear converged at step: " << l << std::endl
          << std::endl;
      d_log(oss);

      break;
    }
  } // nonlinear solver loop

  oss << "\n  End of nonlinear loop\n";
  d_log(oss);
}

void avafv::Model::test_tum_2() {

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

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  oss << "\n  Nonlinear loop\n";
  d_log(oss);

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    oss << "    ____________________\n";
    oss << "    Nonlinear step: " << l << "\n\n";
    oss << "      Solving ";
    d_log(oss);

    // solve nutrient
    oss << "[3D nutrient] -> ";
    d_log(oss);
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    oss << "[tumor species] -> ";
    d_log(oss);
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    oss << "[hypoxic species] -> ";
    d_log(oss);
    hyp.solve();

    // solve necrotic
    oss << "[necrotic species]\n";
    d_log(oss);
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      oss << "      Linear converged at step: " << n_linear_iterations
          << ", residual: " << final_linear_residual
          << ", Nonlinear convergence: ||u - u_old|| = "
          << nonlinear_global_error << std::endl << std::endl;
      d_log(oss);
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      oss << "      Nonlinear converged at step: " << l << std::endl
          << std::endl;
      d_log(oss);

      break;
    }
  } // nonlinear solver loop

  oss << "\n  End of nonlinear loop\n";
  d_log(oss);
}