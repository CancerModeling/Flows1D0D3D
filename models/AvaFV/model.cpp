////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "model.hpp"
#include "systems.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

namespace {
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
                             std::vector<double> &QOI_MASS,
                             const std::string &filename,
                             Parallel::Communicator *comm) {

  // initialize logger
  out << " ********** AvaFV **************\n";

  // init seed for random number
  random_init();

  // read input file
  auto input = avafv::InputDeck(filename);

  // create mesh
  out << "Creating tumor mesh\n";
  ReplicatedMesh mesh(*comm);
  if (input.d_read_mesh_flag)
    mesh.read(input.d_mesh_filename);
  else
    create_mesh(input, mesh);

  out << "Creating tumor system\n";
  EquationSystems tum_sys(mesh);
  // add parameters to system
  tum_sys.parameters.set<avafv::InputDeck *>("input_deck") = &input;
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
  auto model = Model(argc, argv, QOI_MASS, filename, comm, input, mesh, tum_sys,
                     nec, tum, nut, hyp, taf, grad_taf);
}

void avafv::create_mesh(avafv::InputDeck &input, ReplicatedMesh &mesh) {

  double xmax = input.d_domain_params[1];
  double ymax = input.d_domain_params[3];
  double zmax = input.d_domain_params[5];
  unsigned int nelx, nely, nelz;
  double hx, hy, hz;
  if (input.d_dim == 2) {

    if (input.d_use_mesh_size_for_disc) {
      hx = input.d_mesh_size;
      hy = input.d_mesh_size;
      nelx = xmax / hx;
      nely = ymax / hy;
    } else {
      nelx = input.d_num_elems;
      nely = input.d_num_elems;
      hx = xmax  / nelx;
      hy = ymax / nely;
    }

    input.d_elem_face_size = hx;
    input.d_elem_size = hx * hx;
    input.d_face_by_h = 1.;

    // check if length of element in x and y direction are same
    if (std::abs(hx - hy) > 0.001 * hx) {
      libmesh_error_msg("Size of element needs to be same in both direction, "
                        "ie. element needs to be square\n"
                        "If domain is rectangle than specify number of "
                        "elements in x and y different so that element is "
                        "square\n");
    }

    // either read or create mesh
    if (input.d_restart)
      mesh.read(input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_square(mesh, nelx, nely, 0., xmax, 0., ymax,
                                          QUAD4);

  } else if (input.d_dim == 3) {

    if (input.d_use_mesh_size_for_disc) {
      hx = input.d_mesh_size;
      hy = input.d_mesh_size;
      hz = input.d_mesh_size;
      nelx = xmax / hx;
      nely = ymax / hy;
      nelz = zmax / hz;
    } else {
      nelx = input.d_num_elems;
      nely = input.d_num_elems;
      nelz = input.d_num_elems;
      hx = xmax / nelx;
      hy = ymax / nely;
      hz = zmax / nelz;
    }

    input.d_elem_face_size = hx * hx;
    input.d_elem_size = hx * hx * hx;
    input.d_face_by_h = hx;

    // check if length of element in x and y direction are same
    if (std::abs(hx - hy) > 0.001 * hx or std::abs(hx - hz) > 0.001 * hx) {
      libmesh_error_msg("Size of element needs to be same in all three "
                        "direction, ie. element needs to be square\n"
                        "If domain is cuboid than specify number of "
                        "elements in x, y and z different so that element is "
                        "cube\n");
    }

    // either read or create mesh
    if (input.d_restart)
      mesh.read(input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, 0., xmax, 0.,
                                        ymax, 0., zmax, HEX8);
  }
}

// Model class
avafv::Model::Model(
    int argc, char **argv, std::vector<double> &QOI_MASS,
    const std::string &filename, Parallel::Communicator *comm,
    avafv::InputDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
    TransientLinearImplicitSystem &nec, TransientLinearImplicitSystem &tum,
    TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &hyp,
    TransientLinearImplicitSystem &taf, TransientLinearImplicitSystem &grad_taf)
    : d_step(0), d_time(0.), d_dt(0.), d_hmin(0.), d_hmax(0.),
      d_bounding_box(Point(), Point()), d_input(input), d_comm_p(comm),
      d_mesh(mesh), d_tum_sys(tum_sys),
      d_nec_assembly(this, "Necrotic", nec), d_tum_assembly(this, "Tumor", tum),
      d_nut_assembly(this, "Nutrient", nut),
      d_hyp_assembly(this, "Hypoxic", hyp), d_taf_assembly(this, "TAF", taf),
      d_grad_taf_assembly(this, "TAF_Gradient", grad_taf){

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
    // Note: For velocity, we also do not perform assembly
    // vel.attach_assemble_object(d_vel_assembly);

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

    if (d_input.d_perform_output) {
      d_mesh.print_info();
      d_tum_sys.print_info();
      d_mesh.write("mesh.e");
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

  //
  // File to print solution
  //
  // Tumor model
  if (d_input.d_perform_output)
    write_system(0);

  //
  // set time parameters
  //
  d_step = d_input.d_init_step;
  d_time = d_input.d_init_time;
  d_dt = d_input.d_dt;

  d_tum_sys.parameters.set<unsigned int>("linear solver maximum iterations") =
      d_input.d_linear_max_iters;

  //
  // Solve
  //

  do {

    // Prepare time step
    d_step++;
    d_time += d_dt;

    out << std::endl
        << "Solving time step: " << d_step << ", time: " << d_time << std::endl;

    // solve tumor-network system
    // test_tum();
    test_tum_2();

    // Post-processing
    if ((d_step >= d_input.d_dt_output_interval and
         (d_step % d_input.d_dt_output_interval == 0)) or
        d_step == 0) {

      // write tumor solution
      write_system((d_step - d_input.d_init_step) /
                       d_input.d_dt_output_interval,
                   &QOI_MASS);
    }

  } while (d_step < d_input.d_steps);
}

void avafv::Model::write_system(const unsigned int &t_step,
                                 std::vector<double> *QOI_MASS) {

  ExodusII_IO exodus(d_mesh);

  // compute total integral of order parameter u
  double value_mass, total_mass;
  util::computeMass(d_tum_sys, "Tumor", "tumor", value_mass);

  // gather from other processors
  MPI_Allreduce(&value_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  //
  if (QOI_MASS != nullptr)
    QOI_MASS->push_back(total_mass);

  // print to screen
  out << "Time: " << d_time << ", Total tumor Mass: " << total_mass
      << std::endl;

  // print to file
  std::ofstream out_file;
  out_file.precision(16);
  out_file.open("sim_info.txt", std::ios_base::app | std::ios_base::out);
  out_file << std::scientific << d_time << " " << total_mass << std::endl;

  // write mesh and simulation results
  std::string filename = "sim_" + std::to_string(t_step) + ".e";

  // write to exodus
  // exodus.write_timestep(filename, d_tum_sys, 1, d_time);

  //
  VTKIO(d_mesh).write_equation_systems(
      "sim_" + std::to_string(t_step) + ".pvtu", d_tum_sys);

  // save for restart
  if (d_input.d_restart_save &&
      (d_step % d_input.d_dt_restart_save_interval == 0)) {

    out << "Saving files for restart" << std::endl;
    if (t_step == 0) {

      std::string mesh_file = "mesh.e";
      d_mesh.write(mesh_file);
    }

    std::string solutions_file = "solution_" + std::to_string(d_step) + ".e";
    d_tum_sys.write(solutions_file, WRITE);
  }
}

void avafv::Model::project_solution_to_physical_range(
    const MeshBase &mesh, TransientLinearImplicitSystem &sys) {

  if (!d_input.d_project_solution_to_physical_range)
    return;

  // loop over dofs and project solution to [0,1] range
  const auto &sys_name = sys.name();

  unsigned int num_vars = 1;
  unsigned int sys_number = sys.number();
  unsigned int var_number = 0;
  if (sys_name == "Tumor")
    num_vars = 2;

  // loop over nodes and modify dofs
  for (const auto &node : mesh.local_node_ptr_range()) {

    const auto &dof = node->dof_number(sys_number, var_number, 0);

    auto val = sys.current_solution(dof);
    if (val < 0.)
      sys.solution->add(dof, -val);
    else if (val > 1.)
      sys.solution->add(dof, -val + 1.);
    else
      continue;
  }

  sys.solution->close();
  sys.update();
}

void avafv::Model::solve_system() {

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

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;
    out << "Nonlinear step: " << l << "\n";

    // solve nutrient
    out << std::endl << "Solving tissue nutrient system" << std::endl;
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << std::endl << "Solving tumor system" << std::endl;
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << std::endl << "Solving hypoxic system" << std::endl;
    hyp.solve();

    // solve necrotic
    out << std::endl << "Solving necrotic system" << std::endl;
    nec.solve();

    // solve taf
    out << std::endl << "Solving TAF system" << std::endl;
    taf.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "  Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  // solve for gradient of taf
  out << std::endl << "Solving gradient of TAF system" << std::endl;
  //grad_taf.solve();
  d_grad_taf_assembly.solve();
}

void avafv::Model::test_tum() {

  d_test_name = "test_tum";

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

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;
    out << "Nonlinear step: " << l << "\n";

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << std::endl << "Solving tumor system" << std::endl;
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << std::endl << "Solving hypoxic system" << std::endl;
    hyp.solve();

    // solve necrotic
    out << std::endl << "Solving necrotic system" << std::endl;
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "  Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop
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

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;
    out << "Nonlinear step: " << l << "\n";

    // solve nutrient
    out << std::endl << "Solving tissue nutrient system" << std::endl;
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << std::endl << "Solving tumor system" << std::endl;
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << std::endl << "Solving hypoxic system" << std::endl;
    hyp.solve();

    // solve necrotic
    out << std::endl << "Solving necrotic system" << std::endl;
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "  Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop
}