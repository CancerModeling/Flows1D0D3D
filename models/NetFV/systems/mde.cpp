////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfv::initial_condition_mde(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"MDE");

  if (var_name == "mde") {

    const auto *deck = es.get<netfv::InputDeck *>("input_deck");

    return deck->d_mde_ic_val * initial_condition_hyp_kernel(p, deck);
  }

  return 0.;
}

// Assembly class
void netfv::MdeAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfv::MdeAssembly::assemble_face() {

  // call diffusion-advection calculation function
  if (d_model_p->get_input_deck().d_advection_active)
    netfv::assemble_diffusion_advection(d_model_p->get_mde_assembly(),
                             d_model_p->get_pres_assembly(),
                             d_model_p->get_tum_assembly());
  else
    netfv::assemble_diffusion(d_model_p->get_mde_assembly());
}

void netfv::MdeAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real mde_old = 0.;
  Real nut_cur = 0.;
  Real tum_cur = 0.;
  Real nec_cur = 0.;
  Real ecm_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_var_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);

    // const unsigned int n_dofs = mde.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get fields at this element
      nut_cur = nut.get_current_sol(0);
      tum_cur = tum.get_current_sol_var(0, 0);
      nec_cur = nec.get_current_sol(0);
      ecm_cur = ecm.get_current_sol(0);
      mde_old = get_old_sol(0);

      Real aux_1 = deck.d_lambda_MDE_P * (tum_cur - nec_cur) * ecm_cur *
                   deck.d_sigma_HP / (1. + nut_cur);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += mde_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_MDE_D +
                       deck.d_lambda_MDE_P * aux_1 +
                       deck.d_lambda_ECM_D * ecm_cur);
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_MDE_P * aux_1;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void netfv::MdeAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real mde_old = 0.;
  Real nut_proj = 0.;
  Real tum_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_var_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);

    // const unsigned int n_dofs = mde.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get fields at this element
      mde_old = get_old_sol(0);


      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      tum_proj = util::project_concentration(tum.get_current_sol_var(0, 0));
      nec_proj = util::project_concentration(nec.get_current_sol(0));
      ecm_proj = util::project_concentration(ecm.get_current_sol(0));

      Real aux_1 = deck.d_lambda_MDE_P * (tum_proj - nec_proj) * ecm_proj *
                   deck.d_sigma_HP / (1. + nut_proj);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += mde_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_MDE_D +
                        aux_1 +
                       deck.d_lambda_ECM_D * ecm_proj);
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * aux_1;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void netfv::MdeAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real mde_old = 0.;
  Real mde_cur = 0.;
  Real nut_proj = 0.;
  Real tum_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_var_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);

    // const unsigned int n_dofs = mde.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get fields at this element
      mde_old = get_old_sol(0);
      mde_cur = get_current_sol(0);

      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      tum_proj = util::project_concentration(tum.get_current_sol_var(0, 0));
      nec_proj = util::project_concentration(nec.get_current_sol(0));
      ecm_proj = util::project_concentration(ecm.get_current_sol(0));

      Real aux_1 = deck.d_lambda_MDE_P * (tum_proj - nec_proj) * ecm_proj *
                   deck.d_sigma_HP / (1. + nut_proj);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += mde_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_MDE_D +
                       deck.d_lambda_MDE_P * aux_1 +
                       deck.d_lambda_ECM_D * ecm_proj);
      Fe(0) += -a_source * mde_cur;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_MDE_P * aux_1;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}