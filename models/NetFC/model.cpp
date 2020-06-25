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
#include "rw/vtk_io.hpp"
#include <random>

namespace {
void random_init() { srand(time(nullptr)); }
} // namespace

namespace netfc {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  double init_time = es.parameters.get<Real>("init_time");
  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);
  es.parameters.set<Real>("time") = sys.time = init_time;

  if (system_name == "Nutrient")
    sys.project_solution(netfc::initial_condition_nut, nullptr, es.parameters);
  else if (system_name == "Tumor")
    sys.project_solution(netfc::initial_condition_tum, nullptr, es.parameters);
  else if (system_name == "Hypoxic")
    sys.project_solution(netfc::initial_condition_hyp, nullptr, es.parameters);
  else if (system_name == "Necrotic")
    sys.project_solution(netfc::initial_condition_nec, nullptr, es.parameters);
  else if (system_name == "TAF")
    sys.project_solution(netfc::initial_condition_taf, nullptr, es.parameters);
  else if (system_name == "ECM")
    sys.project_solution(netfc::initial_condition_ecm, nullptr, es.parameters);
  else if (system_name == "MDE")
    sys.project_solution(netfc::initial_condition_mde, nullptr, es.parameters);
  else if (system_name == "Pressure")
    sys.project_solution(netfc::initial_condition_pres, nullptr,
                         es.parameters);
  else {
    return;
  }
}
} // namespace netfc

// Model class
netfc::Model::Model(int argc, char **argv, std::vector<double> &QOI_MASS,
                     const std::string &filename, Parallel::Communicator *comm)
    : d_step(0), d_time(0.), d_dt(0.), d_hmin(0.), d_hmax(0.),
      d_bounding_box(Point(), Point()), d_comm_p(comm),
      d_mesh(ReplicatedMesh(*d_comm_p)), d_tum_sys(d_mesh), d_network(this),
      d_nec_assembly(this, "Necrotic"), d_tum_assembly(this, "Tumor"),
      d_nut_assembly(this, "Nutrient"), d_hyp_assembly(this, "Hypoxic"),
      d_taf_assembly(this, "TAF"), d_ecm_assembly(this, "ECM"),
      d_mde_assembly(this, "MDE"), d_pres_assembly(this, "Pressure"),
      d_grad_taf_assembly(this, "TAF_Gradient"),
      d_vel_assembly(this, "Velocity") {

  // initialize logger
  out << " ********** NetTum **************\n";

  // init seed for random number
  random_init();

  // read input file
  d_input = netfc::InputDeck(filename);
  // d_input.print();
  // bounding box
  d_bounding_box.first =
      Point(d_input.d_domain_params[0], d_input.d_domain_params[2],
            d_input.d_domain_params[4]);
  d_bounding_box.second =
      Point(d_input.d_domain_params[1], d_input.d_domain_params[3],
            d_input.d_domain_params[5]);

  //
  // 2d/3d tumor growth model
  //
  // create mesh
  out << "Creating tumor mesh\n";
  if (d_input.d_read_mesh_flag)
    d_mesh.read(d_input.d_mesh_filename);
  else
    create_mesh();

  // get point locator
  const auto &mesh_locator = d_mesh.point_locator();

  out << "Creating 1-D network\n";
  d_network.create_initial_network();

  out << "Creating tumor system\n";
  setup_system();

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

  do{

    // Prepare time step
    d_step++;
    d_time += d_dt;

    out << " " << std::endl;
    out << "-----------------------------------------------------" << std::endl;
    out << "Solving time step: " << d_step << ", time: " << d_time << std::endl;
    out << "-----------------------------------------------------" << std::endl;

    d_network.solve3D1DFlowProblem( d_step, d_time );
    d_network.solve3D1DNutrientProblem( d_step, d_time );
   // d_network.solve3DProlificCellProblem( d_step, d_time );
    d_network.solve3DTAFProblem( d_step, d_time );

    if( d_step%2 ){

        d_network.updateNetwork();

    }

    //d_network.printDataVGM();

  }while (d_step < d_input.d_steps);

}

void netfc::Model::create_mesh() {

  if (d_input.d_dim == 2) {

    unsigned int nelx = d_input.d_num_elems_vec[0];
    unsigned int nely = d_input.d_num_elems_vec[1];
    double xmax = d_input.d_domain_params[1];
    double ymax = d_input.d_domain_params[3];

    d_input.d_mesh_size_vec[0] = (xmax - 0.) / nelx;
    d_input.d_mesh_size_vec[1] = (ymax - 0.) / nely;

    // check if length of element in x and y direction are same
    if (std::abs(d_input.d_mesh_size_vec[0] - d_input.d_mesh_size_vec[1]) >
        0.001 * d_input.d_mesh_size_vec[0]) {
      libmesh_error_msg("Size of element needs to be same in both direction, "
                        "ie. element needs to be square\n"
                        "If domain is rectangle than specify number of "
                        "elements in x and y different so that element is "
                        "square\n");
    }

    // either read or create mesh
    if (d_input.d_restart)
      d_mesh.read(d_input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_square(d_mesh, nelx, nely, 0., xmax, 0.,
                                          ymax, QUAD4);

  } else if (d_input.d_dim == 3) {

    unsigned int nelx = d_input.d_num_elems_vec[0];
    unsigned int nely = d_input.d_num_elems_vec[1];
    unsigned int nelz = d_input.d_num_elems_vec[2];
    double xmax = d_input.d_domain_params[1];
    double ymax = d_input.d_domain_params[3];
    double zmax = d_input.d_domain_params[5];

    d_input.d_mesh_size_vec[0] = (xmax - 0.) / nelx;
    d_input.d_mesh_size_vec[1] = (ymax - 0.) / nely;
    d_input.d_mesh_size_vec[2] = (zmax - 0.) / nelz;

    // check if length of element in x and y direction are same
    if (std::abs(d_input.d_mesh_size_vec[0] - d_input.d_mesh_size_vec[1]) >
        0.001 * d_input.d_mesh_size_vec[0] or
        std::abs(d_input.d_mesh_size_vec[0] - d_input.d_mesh_size_vec[2]) >
        0.001 * d_input.d_mesh_size_vec[0]) {
      libmesh_error_msg("Size of element needs to be same in all three "
                        "direction, ie. element needs to be square\n"
                        "If domain is cuboid than specify number of "
                        "elements in x, y and z different so that element is "
                        "cube\n");
    }

    // either read or create mesh
    if (d_input.d_restart)
      d_mesh.read(d_input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_cube(d_mesh, nelx, nely, nelz, 0., xmax, 0.,
                                        ymax, 0., zmax, HEX8);
  }

  // modify some values in input deck
  d_hmax = d_input.d_mesh_size_vec[0];
  d_hmin = d_input.d_mesh_size_vec[0];
}

void netfc::Model::setup_system() {

  // add parameters to system
  d_tum_sys.parameters.set<netfc::InputDeck *>("input_deck") = &d_input;
  d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = 0;
  d_tum_sys.parameters.set<Real>("init_time") = d_input.d_init_time;
  d_tum_sys.parameters.set<Real>("time_step") = d_input.d_dt;

  // read if available
  if (d_input.d_restart)
    d_tum_sys.read(d_input.d_sol_restart_file, READ);

  // Add systems, variables and assemble
  auto &nut = d_tum_sys.add_system<TransientLinearImplicitSystem>("Nutrient");
  auto &tum = d_tum_sys.add_system<TransientLinearImplicitSystem>("Tumor");
  auto &hyp = d_tum_sys.add_system<TransientLinearImplicitSystem>("Hypoxic");
  auto &nec = d_tum_sys.add_system<TransientLinearImplicitSystem>("Necrotic");
  auto &taf = d_tum_sys.add_system<TransientLinearImplicitSystem>("TAF");
  auto &ecm = d_tum_sys.add_system<TransientLinearImplicitSystem>("ECM");
  auto &mde = d_tum_sys.add_system<TransientLinearImplicitSystem>("MDE");
  auto &pres = d_tum_sys.add_system<TransientLinearImplicitSystem>("Pressure");
  auto &grad_taf =
      d_tum_sys.add_system<TransientLinearImplicitSystem>("TAF_Gradient");
  auto &vel = d_tum_sys.add_system<TransientLinearImplicitSystem>("Velocity");

  if (d_input.d_restart) {
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

  if (!d_input.d_restart) {
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
    if (d_input.d_dim > 2)
      grad_taf.add_variable("taf_gradz", FIRST);

    // variable in gradient TAF system
    vel.add_variable("velocity_x", FIRST);
    vel.add_variable("velocity_y", FIRST);
    if (d_input.d_dim > 2)
      vel.add_variable("velocity_z", FIRST);

    // attach initial condition function to systems
    tum.attach_init_function(initial_condition);
    hyp.attach_init_function(initial_condition);
    nut.attach_init_function(initial_condition);
    ecm.attach_init_function(initial_condition);
    mde.attach_init_function(initial_condition);
    pres.attach_init_function(initial_condition);
  }

  //
  // attach assembly function
  //
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

  // Add boundary condition
  boundary_condition_nut(d_tum_sys);

  // For pressure, we explicitly constrain the matrix
  // boundary_condition_pres(d_tum_sys);

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

}

void netfc::Model::write_system(const unsigned int &t_step,
                                 std::vector<double> *QOI_MASS) {

  ExodusII_IO exodus(d_mesh);

  auto &pres = d_tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");
  std::vector<Number> p_save;
  std::vector<unsigned int> p_dofs;
  {
    double factor = d_input.d_mmhgFactor;
    // double factor = 1.0;

    const unsigned int v_pres = pres.variable_number("pressure");
    const DofMap &pres_map = pres.get_dof_map();
    std::vector<unsigned int> dof_indices_pres;

    // Looping through elements
    MeshBase::const_element_iterator el = d_mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
        d_mesh.active_local_elements_end();

    for (; el != end_el; ++el) {

      const Elem *elem = *el;
      pres_map.dof_indices(elem, dof_indices_pres, v_pres);

      double pt = pres.current_solution(dof_indices_pres[0]);

      p_save.push_back(pt);
      p_dofs.push_back(dof_indices_pres[0]);

      pt = pt / factor;

      pres.solution->set(dof_indices_pres[0], pt);
    }

    pres.solution->close();
    pres.update();
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
  out << "Time: " << d_time << ", Total tumor Mass: " << total_mass
      << std::endl;

  // print to file
  std::ofstream out_file;
  out_file.precision(16);
  out_file.open("sim_info.txt", std::ios_base::app | std::ios_base::out);
  out_file << std::scientific << d_time << " " << total_mass << std::endl;

  // write mesh and simulation results
  std::string filename = "sim_" +std::to_string(t_step) + ".e";

  // write to exodus
  // exodus.write_timestep(filename, d_tum_sys, 1, d_time);

  //
  rw::VTKIO(d_mesh).write_equation_systems("sim_" + std::to_string(t_step) +
                                       ".pvtu",  d_tum_sys);

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

  {
    for (unsigned int i = 0; i < p_dofs.size(); i++) {

      pres.solution->set(p_dofs[i], p_save[i]);
    }
    pres.solution->close();
    pres.update();
  }
}

void netfc::Model::project_solution_to_physical_range(
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

void netfc::Model::solve_system() {

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
  auto &vel = d_tum_sys.get_system<TransientLinearImplicitSystem>("Velocity");

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
  *pres.old_local_solution = *pres.current_local_solution;
  *vel.old_local_solution = *vel.current_local_solution;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_tum(
      tum.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;

    // solver for 1-D pressure and nutrient
    out << std::endl << "Solving Network system" << std::endl;
    //d_network.solve_system();

    // solve for pressure in tissue
    out << std::endl << "Solving tissue pressure system" << std::endl;
    pres.solve();

    // solve for velocity in tissue
    out << std::endl << "Solving tissue velocity system" << std::endl;
    vel.solve();

    // solve nutrient<
    out << std::endl << "Solving nutrient system" << std::endl;
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

    // solve mde
    out << std::endl << "Solving MDE system" << std::endl;
    mde.solve();

    // solve ecm
    out << std::endl << "Solving ECM system" << std::endl;
    ecm.solve();

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
  grad_taf.solve();
}

void netfc::Model::test_nutrient() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

  // update time
  nut.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_nut(
      nut.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // Nutrient coupling
  out << "Nonlinear loop for coupled nutrient systems\n";

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
    //for (unsigned int l = 0; l < 1; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;

    // solver for 1-D nutrient
    out << std::endl << "Solving network nutrient system" << std::endl;
    //d_network.solveVGMforNutrient();

    // solve for nutrient in tissue
    last_nonlinear_soln_nut->zero();
    last_nonlinear_soln_nut->add(*nut.solution);
    out << std::endl << "Solving tissue nutrient system" << std::endl;
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

      std::cout << "  Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error_nut << std::endl;
    }
    if (nonlinear_global_error_nut < d_input.d_nonlin_tol) {

      std::cout << "Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }

  } // nonlinear solver loop
}

void netfc::Model::solve_nutrient_3D(){

     out << " " << std::endl;
     out << "Solve 3D nutrient equation" << std::endl;
     
     // get systems
     auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

     // update time
     nut.time = d_time;
     d_tum_sys.parameters.set<Real>("time") = d_time;

     // update old solution
     *nut.old_local_solution = *nut.current_local_solution;
     nut.solve();

}

void netfc::Model::test_nutrient_2() {

  // get systems
  auto &nut = d_tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");

  // update time
  nut.time = d_time;
  d_tum_sys.parameters.set<Real>("time") = d_time;

  // update old solution
  *nut.old_local_solution = *nut.current_local_solution;

  // solver for 1-D nutrient
  out << std::endl << "Solving network nutrient system" << std::endl;
  //d_network.solveVGMforNutrient();

  out << std::endl << "Solving tissue nutrient system" << std::endl;
  nut.solve();

  return;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln_nut(
      nut.solution->clone());

  // Nonlinear iteration loop
  d_tum_sys.parameters.set<Real>("linear solver tolerance") =
      d_input.d_linear_tol;

  // Nutrient coupling
  out << "Nonlinear loop for 3d nutrient systems\n";

  // nonlinear loop
  for (unsigned int l = 0; l < d_input.d_nonlin_max_iters; ++l) {
    //for (unsigned int l = 0; l < 1; ++l) {

    d_tum_sys.parameters.set<unsigned int>("nonlinear_step") = l;

    // solver for 1-D nutrient
    out << std::endl << "Solving network nutrient system" << std::endl;
    //d_network.solveVGMforNutrient();

    // solve for nutrient in tissue
    last_nonlinear_soln_nut->zero();
    last_nonlinear_soln_nut->add(*nut.solution);
    out << std::endl << "Solving tissue nutrient system" << std::endl;
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

      std::cout << "  Linear converged at step: " << n_linear_iterations
                << ", residual: " << final_linear_residual
                << ", Nonlinear convergence: ||u - u_old|| = "
                << nonlinear_global_error_nut << std::endl;
    }
    if (nonlinear_global_error_nut < d_input.d_nonlin_tol) {

      std::cout << "Nonlinear converged at step: " << l << std::endl
                << std::endl;

      break;
    }

  } // nonlinear solver loop
}

void netfc::Model::test_taf() {

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

  // update time
  nut.time = d_time;
  tum.time = d_time;
  hyp.time = d_time;
  nec.time = d_time;
  taf.time = d_time;
  ecm.time = d_time;
  mde.time = d_time;
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

    out << std::endl << "Solving network nutrient system" << std::endl;
    //d_network.solveVGMforNutrient();

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

    // solve mde
    out << std::endl << "Solving MDE system" << std::endl;
    // mde.solve();

    // solve ecm
    out << std::endl << "Solving ECM system" << std::endl;
    // ecm.solve();

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
  grad_taf.solve();
}

void netfc::Model::test_taf_2() {

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

  out << std::endl << "Solving TAF system" << std::endl;
  taf.solve();

  // solve for gradient of taf
  out << std::endl << "Solving gradient of TAF system" << std::endl;
  grad_taf.solve();
}
