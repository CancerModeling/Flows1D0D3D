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
#include "ghosting_functor.hpp"
#include <random>

namespace {
void random_init() { srand(time(nullptr)); }
} // namespace

namespace netfv {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  if (system_name == "Nutrient")
    sys.project_solution(netfv::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Tumor")
    sys.project_solution(netfv::initial_condition_tum, nullptr, es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(netfv::initial_condition_hyp, nullptr, es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(netfv::initial_condition_nec, nullptr, es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(netfv::initial_condition_taf, nullptr, es.parameters);
  else if (system_name == "ECM")
    sys.project_solution(netfv::initial_condition_ecm, nullptr, es.parameters);
  else if (system_name == "MDE")
    sys.project_solution(netfv::initial_condition_mde, nullptr, es.parameters);
  else {
    return;
  }
}
} // namespace netfv

// Model setup and run
void netfv::model_setup_run(int argc, char **argv,
                             std::vector<double> &QOI_MASS,
                             const std::string &filename,
                             Parallel::Communicator *comm) {

  // initialize logger
  out << " ********** NetTum **************\n";

  // init seed for random number
  random_init();

  // read input file
  auto input = netfv::InputDeck(filename);

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
  tum_sys.parameters.set<netfv::InputDeck *>("input_deck") = &input;
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
      tum.add_variable("tumor", CONSTANT, MONOMIAL);
      tum.add_variable("chemical_tumor", CONSTANT, MONOMIAL);

      // variable in hypoxic system
      hyp.add_variable("hypoxic", CONSTANT, MONOMIAL);

      // variable in necrotic system
      nec.add_variable("necrotic", CONSTANT, MONOMIAL);

      // variable in TAF system
      taf.add_variable("taf", CONSTANT, MONOMIAL);

      // variable in ECM system
      ecm.add_variable("ecm", CONSTANT, MONOMIAL);

      // variable in MDE system
      mde.add_variable("mde", CONSTANT, MONOMIAL);

      // variable in Pressure system
      pres.add_variable("pressure", CONSTANT, MONOMIAL);

      // variable in gradient TAF system
      grad_taf.add_variable("taf_gradx", CONSTANT, MONOMIAL);
      grad_taf.add_variable("taf_grady", CONSTANT, MONOMIAL);
      if (input.d_dim > 2)
        grad_taf.add_variable("taf_gradz", CONSTANT, MONOMIAL);

      // variable in velocity system
      vel.add_variable("velocity_x", CONSTANT, MONOMIAL);
      vel.add_variable("velocity_y", CONSTANT, MONOMIAL);
      if (input.d_dim > 2)
        vel.add_variable("velocity_z", CONSTANT, MONOMIAL);

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

  // create ghosting functor
  //  netfv::GhostingFunctorNet ghosting_fun(mesh);
  //  tum.get_dof_map().add_coupling_functor(ghosting_fun);
  //  hyp.get_dof_map().add_coupling_functor(ghosting_fun);
  //  nec.get_dof_map().add_coupling_functor(ghosting_fun);
  //  taf.get_dof_map().add_coupling_functor(ghosting_fun);
  //  mde.get_dof_map().add_coupling_functor(ghosting_fun);
  //  ecm.get_dof_map().add_coupling_functor(ghosting_fun);

  //
  // Create Model class
  //
  auto model = Model(argc, argv, QOI_MASS, filename, comm, input, mesh, tum_sys,
                     nec, tum, nut, hyp, taf, ecm, mde, pres, grad_taf, vel);
}

void netfv::create_mesh(netfv::InputDeck &input, ReplicatedMesh &mesh) {

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
netfv::Model::Model(
    int argc, char **argv, std::vector<double> &QOI_MASS,
    const std::string &filename, Parallel::Communicator *comm,
    netfv::InputDeck &input, ReplicatedMesh &mesh, EquationSystems &tum_sys,
    TransientLinearImplicitSystem &nec, TransientLinearImplicitSystem &tum,
    TransientLinearImplicitSystem &nut, TransientLinearImplicitSystem &hyp,
    TransientLinearImplicitSystem &taf, TransientLinearImplicitSystem &ecm,
    TransientLinearImplicitSystem &mde, TransientLinearImplicitSystem &pres,
    TransientLinearImplicitSystem &grad_taf,
    TransientLinearImplicitSystem &vel)
    : d_step(0), d_time(0.), d_dt(0.), d_hmin(0.), d_hmax(0.),
      d_bounding_box(Point(), Point()), d_nonlinear_step(0),
      d_is_output_step(false),
      d_input(input), d_comm_p(comm),
      d_mesh(mesh), d_tum_sys(tum_sys), d_network(this),
      d_nec_assembly(this, "Necrotic", nec), d_tum_assembly(this, "Tumor", tum),
      d_nut_assembly(this, "Nutrient", nut),
      d_hyp_assembly(this, "Hypoxic", hyp), d_taf_assembly(this, "TAF", taf),
      d_ecm_assembly(this, "ECM", ecm), d_mde_assembly(this, "MDE", mde),
      d_pres_assembly(this, "Pressure", pres),
      d_grad_taf_assembly(this, "TAF_Gradient", grad_taf),
      d_vel_assembly(this, "Velocity", vel){

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
    ecm.time = d_input.d_init_time;
    mde.time = d_input.d_init_time;
    pres.time = d_input.d_init_time;
    grad_taf.time = d_input.d_init_time;
    vel.time = d_input.d_init_time;

    if (d_input.d_perform_output) {
      d_mesh.print_info();
      d_tum_sys.print_info();
      d_mesh.write("mesh.e");
    }

    // set Petsc matrix option to suppress the error
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(pres.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
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
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(mde.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    {
      PetscMatrix<Number> *pet_mat =
          dynamic_cast<PetscMatrix<Number> *>(ecm.matrix);
      MatSetOption(pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
  }

  // 1-D network
  out << "\n\nCreating 1-D network\n";
  d_network.create_initial_network();

  //
  // File to print solution
  //
  // Tumor model
  if (d_input.d_perform_output)
    write_system(0);

  // network
  d_network.writeDataToVTKTimeStep_VGM(0);

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

  //  if (d_input.d_test_name == "test_net_tum_2") {
  //    // decouple 1d-3d pressure
  //    d_input.d_tissue_flow_L_p = 0.;
  //
  //    // solve only for 1d pressure
  //    out << std::endl << "Solving Network pressure system" << std::endl;
  //    d_network.solveVGMforPressure();
  //  }

  if (!d_input.d_test_name.empty())
    out << "\n\nSolving sub-system for test: " << d_input.d_test_name << "\n\n";

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
    write_system(1, &QOI_MASS);
    d_network.writeDataToVTKTimeStep_VGM(1);
    return;
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

    out << "\n\n________________________________________________________\n";
    out << "At time step: " << d_step << ", time: " << d_time << "\n\n";

    // solve tumor-network system
    // out << std::endl << "Solving Tumor-Network coupled system" << std::endl;
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
      solve_system();

    // update network
    if (d_step % d_input.d_network_update_interval == 0) {
      out << "\n  Updating Network\n";
      d_network.update_network();
    }

    // Post-processing
    if (d_is_output_step) {

      // write tumor solution
      write_system((d_step - d_input.d_init_step) /
                       d_input.d_dt_output_interval,
                   &QOI_MASS);
      d_network.writeDataToVTKTimeStep_VGM((d_step - d_input.d_init_step) /
                                           d_input.d_dt_output_interval);
    }

  } while (d_step < d_input.d_steps);
}

void netfv::Model::write_system(const unsigned int &t_step,
                                 std::vector<double> *QOI_MASS) {

  ExodusII_IO exodus(d_mesh);

  std::vector<Number> p_save;
  std::vector<unsigned int> p_dofs;
  {
    double factor = d_input.d_mmhgFactor;
    // double factor = 1.0;

    // Looping through elements
    for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

      d_pres_assembly.init_dof(elem);
      double pt = d_pres_assembly.get_current_sol(0);

      p_save.push_back(pt);
      p_dofs.push_back(d_pres_assembly.get_global_dof_id(0));

      pt = pt / factor;

      d_pres_assembly.d_sys.solution->set(d_pres_assembly.get_global_dof_id(0),
                                          pt);
    }

    d_pres_assembly.d_sys.solution->close();
    d_pres_assembly.d_sys.update();
  }

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
  out << "\n\n  Total tumor Mass: " << total_mass << "\n";

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

    out << "  Saving files for restart" << std::endl;
    if (t_step == 0) {

      std::string mesh_file = "mesh.e";
      d_mesh.write(mesh_file);
    }

    std::string solutions_file = "solution_" + std::to_string(d_step) + ".e";
    d_tum_sys.write(solutions_file, WRITE);
  }

  {
    for (unsigned int i = 0; i < p_dofs.size(); i++) {

      d_pres_assembly.d_sys.solution->set(p_dofs[i], p_save[i]);
    }
    d_pres_assembly.d_sys.solution->close();
    d_pres_assembly.d_sys.update();
  }
}

void netfv::Model::project_solution_to_physical_range(
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

void netfv::Model::solve_system() {

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

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  out << "\n  Nonlinear loop\n\n";

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solver for 1-D pressure and nutrient
    out << "[1D nutrient] -> ";
    d_network.solveVGMforNutrient();

    // solve nutrient
    out << "[3D nutrient] -> ";
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species] -> ";
    nec.solve();

    // solve taf
    out << "[taf species] -> ";
    taf.solve();

    // solve mde
    out << "[mde species] -> ";
    mde.solve();

    // solve ecm
    out << "[ecm species]\n";
    ecm.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";

  // solve for gradient of taf
  out << "      Solving [gradient of taf] -> ";
  d_grad_taf_assembly.solve();

  // solve for velocity
  out << "[velocity]\n";
  d_vel_assembly.solve();
}

void netfv::Model::test_nut() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

  // update time
  nut.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_nut(
      nut.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  out << "\n  Nonlinear loop\n";

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
    // for (unsigned int l = 0; l < 1; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solver for 1-D nutrient
    out << "[1D nutrient] -> ";
    d_network.solveVGMforNutrient();

    // solve for nutrient in tissue
    last_nonlinear_soln_nut->zero();
    last_nonlinear_soln_nut->add(*nut.solution);
    out << "[3D nutrient]\n";
    nut.solve();
    last_nonlinear_soln_nut->add(-1., *nut.solution);
    last_nonlinear_soln_nut->close();

    double nonlinear_loc_error_nut = last_nonlinear_soln_nut->linfty_norm();
    double nonlinear_global_error_nut = 0.;
    MPI_Allreduce(&nonlinear_loc_error_nut, &nonlinear_global_error_nut, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      // const unsigned int n_linear_iterations = tum.n_linear_iterations();
      // const Real final_linear_residual = tum.final_linear_residual();
      const unsigned int n_linear_iterations = nut.n_linear_iterations();
      const Real final_linear_residual = nut.final_linear_residual();

      std::cout << "    Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error_nut << std::endl << std::endl;
    }
    if (nonlinear_global_error_nut < d_input.d_nonlin_tol) {

      std::cout << "    Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }

  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";
}

void netfv::Model::test_nut_2() {

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
  out << "      Solving [1D nutrient] -> ";
  d_network.solveVGMforNutrient();

  out << "[3D nutrient]\n";
  nut.solve();
}

void netfv::Model::solve_pressure() {
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

  // Solve pressure system
  out << "\n\n  Nonlinear loop for pressure\n";
  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
  // Debug
  // for (unsigned int l = 0; l < 10; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solver for 1-D pressure and nutrient
    out << "[1D pressure] -> ";
    d_network.solveVGMforPressure();

    // solve for pressure in tissue
    last_nonlinear_soln_pres->zero();
    last_nonlinear_soln_pres->add(*pres.solution);
    out << "[3D pressure]\n";
    pres.solve();
    last_nonlinear_soln_pres->add(-1., *pres.solution);
    last_nonlinear_soln_pres->close();

    // Nonlinear iteration error
    double nonlinear_loc_error_pres = last_nonlinear_soln_pres->linfty_norm();
    double nonlinear_global_error_pres = 0.;
    MPI_Allreduce(&nonlinear_loc_error_pres, &nonlinear_global_error_pres, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = pres.n_linear_iterations();
      const Real final_linear_residual = pres.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error_pres << std::endl << std::endl;
    }
    if (nonlinear_global_error_pres < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }

  } // nonlinear solver loop

  out << "\n  End of nonlinear loop for pressure\n";
}

void netfv::Model::test_taf() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &ecm = d_tum_sys.get_system<TransientLinearImplicitSystem>("ECM");
  auto &mde = d_tum_sys.get_system<TransientLinearImplicitSystem>("MDE");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");
  auto &vel =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("Velocity");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  ecm.time = d_time;
  mde.time = d_time;
  grad_taf.time = d_time;
  vel.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;
  *taf.old_local_solution = *taf.current_local_solution;
  *ecm.old_local_solution = *ecm.current_local_solution;
  *mde.old_local_solution = *mde.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // solve for pressure
  solve_pressure();

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  out << "\n  Nonlinear loop\n";

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    out << "[1D nutrient] -> ";
    d_network.solveVGMforNutrient();

    out << "[3D nutrient] -> ";
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species] -> ";
    nec.solve();

    // solve taf
    out << "[taf species] -> ";
    taf.solve();

    // solve mde
    out << "[mde species] -> ";
    // mde.solve();

    // solve ecm
    out << "[ecm species]\n";
    // ecm.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";

  // solve for gradient of taf
  out << "      Solving [gradient of taf] -> ";
  //grad_taf.solve();
  d_grad_taf_assembly.solve();

  // solve for gradient of taf
  out << "[velocity]\n";
  d_vel_assembly.solve();
}

void netfv::Model::test_taf_2() {

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

  out << "      Solving [taf] -> ";
  taf.solve();

  // solve for gradient of taf
  out << "[gradient of taf]\n";
  d_grad_taf_assembly.solve();
}

void netfv::Model::test_tum() {

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

  out << "\n  Nonlinear loop\n";

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species]\n";
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";
}

void netfv::Model::test_tum_2() {

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

  out << "\n  Nonlinear loop\n";

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solve nutrient
    out << "[3D nutrient] -> ";
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species]\n";
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";
}

void netfv::Model::test_net_tum() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &vel =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("Velocity");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  vel.time = d_time;
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

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  out << "\n  Nonlinear loop\n";

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solver for 1-D pressure and nutrient
    out << "[1D nutrient] -> ";
    d_network.solveVGMforNutrient();

    // solve nutrient
    out << "[3D nutrient] -> ";
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species]\n";
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";

  // compute below only when we are performing output as these do not play
  // direct role in evolution of sub-system
  if (d_is_output_step) {

    // solve for velocity
    out << "      Solving [velocity] -> ";
    d_vel_assembly.solve();

    // solve for taf
    *taf.old_local_solution = *taf.current_local_solution;
    out << "[taf species] -> ";
    taf.solve();

    // solve for gradient of taf
    out << "[gradient of taf]\n";
    d_grad_taf_assembly.solve();
  }
}

void netfv::Model::test_net_tum_2() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.get_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.get_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.get_system<TransientLinearImplicitSystem>("Necrotic");
  auto &vel =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("Velocity");
  auto &taf = d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF");
  auto &grad_taf =
      d_tum_sys.get_system<TransientLinearImplicitSystem>("TAF_Gradient");

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  vel.time = d_time;
  grad_taf.time = d_time;

  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;
  *tum.old_local_solution = *tum.current_local_solution;
  *hyp.old_local_solution = *hyp.current_local_solution;
  *nec.old_local_solution = *nec.current_local_solution;

  // update old concentration in network
  d_network.update_old_concentration();

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  out << "\n  Nonlinear loop\n";

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_nonlinear_step = l;
    out << "    ____________________\n";
    out << "    Nonlinear step: " << l << "\n\n";
    out << "      Solving ";

    // solver for 1-D pressure and nutrient
    out << "[1D nutrient] -> ";
    d_network.solveVGMforNutrient();

    // solve nutrient
    out << "[3D nutrient] -> ";
    nut.solve();

    // solve tumor
    last_nonlinear_soln_tum->zero();
    last_nonlinear_soln_tum->add(*tum.solution);
    out << "[tumor species] -> ";
    tum.solve();
    last_nonlinear_soln_tum->add(-1., *tum.solution);
    last_nonlinear_soln_tum->close();

    // solve hypoxic
    out << "[hypoxic species] -> ";
    hyp.solve();

    // solve necrotic
    out << "[necrotic species]\n";
    nec.solve();

    // Nonlinear iteration error
    double nonlinear_loc_error = last_nonlinear_soln_tum->linfty_norm();
    double nonlinear_global_error = 0.;
    MPI_Allreduce(&nonlinear_loc_error, &nonlinear_global_error, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    if (d_input.d_perform_output) {

      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();

      std::cout << "      Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error << std::endl << std::endl;
    }
    if (nonlinear_global_error < d_input.d_nonlin_tol) {

      std::cout << "      Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }
  } // nonlinear solver loop

  out << "\n  End of nonlinear loop\n";

  // compute below only when we are performing output as these do not play
  // direct role in evolution of sub-system
  if (d_is_output_step) {

    // solve for velocity
    out << "      Solving [velocity] -> ";
    d_vel_assembly.solve();

    // solve for taf
    *taf.old_local_solution = *taf.current_local_solution;
    out << "[taf species] -> ";
    taf.solve();

    // solve for gradient of taf
    out << "[gradient of taf]\n";
    d_grad_taf_assembly.solve();
  }
}

