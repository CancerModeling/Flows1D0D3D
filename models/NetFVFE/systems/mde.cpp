////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "systems.hpp"
#include "../model.hpp"

Number netfvfe::initial_condition_mde(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"MDE");

  if (var_name == "mde") {

    const auto *deck = es.get<netfvfe::InputDeck *>("input_deck");

    return deck->d_mde_ic_val * initial_condition_hyp_kernel(p, deck);
  }

  return 0.;
}

void netfvfe::boundary_condition_mde(EquationSystems &es) {

  const auto *deck = es.parameters.get<netfvfe::InputDeck *>("input_deck");
}

void netfvfe::MdeAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfvfe::MdeAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // ECM system
  auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
  const unsigned int v_ecm = ecm.variable_number("ecm");
  const DofMap &ecm_map = ecm.get_dof_map();
  std::vector<unsigned int> dof_indices_ecm;

  // MDE system
  auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
  const unsigned int v_mde = mde.variable_number("mde");
  const DofMap &mde_map = mde.get_dof_map();
  std::vector<unsigned int> dof_indices_mde;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = mde.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfvfe::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_mde.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number mde_cur = 0.;
      Number mde_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Number ecm_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
        mde_old += phi[l][qp] * mde.old_solution(dof_indices_mde[l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
      }

      Number aux_1 = deck->d_lambda_MDE_P * (tum_cur - nec_cur) *
                     ecm_cur * deck->d_sigma_HP / (1. + nut_cur);

      Number compute_rhs = JxW[qp] * (mde_old + dt * aux_1);

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_MDE_D +
                     dt * aux_1 + dt * deck->d_lambda_ECM_D * ecm_cur);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_MDE * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    mde_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_mde);
    mde.matrix->add_matrix(Ke, dof_indices_mde);
    mde.rhs->add_vector(Fi, dof_indices_mde);
  }
}

void netfvfe::MdeAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // ECM system
  auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
  const unsigned int v_ecm = ecm.variable_number("ecm");
  const DofMap &ecm_map = ecm.get_dof_map();
  std::vector<unsigned int> dof_indices_ecm;

  // MDE system
  auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
  const unsigned int v_mde = mde.variable_number("mde");
  const DofMap &mde_map = mde.get_dof_map();
  std::vector<unsigned int> dof_indices_mde;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = mde.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfvfe::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_mde.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number mde_cur = 0.;
      Number mde_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Number ecm_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
        mde_old += phi[l][qp] * mde.old_solution(dof_indices_mde[l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
      }

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number nec_proj = util::project_concentration(nec_cur);
      Number ecm_proj = util::project_concentration(ecm_cur);
      Number mde_proj = util::project_concentration(mde_cur);

      Number aux_1 = deck->d_lambda_MDE_P * (tum_proj - nec_proj) *
                     ecm_proj * deck->d_sigma_HP / (1. + nut_proj);

      Number compute_rhs = JxW[qp] * (mde_old + dt * aux_1);

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_MDE_D +
                     dt * aux_1 + dt * deck->d_lambda_ECM_D * ecm_proj);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_MDE * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    mde_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_mde);
    mde.matrix->add_matrix(Ke, dof_indices_mde);
    mde.rhs->add_vector(Fi, dof_indices_mde);
  }
}

void netfvfe::MdeAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // ECM system
  auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
  const unsigned int v_ecm = ecm.variable_number("ecm");
  const DofMap &ecm_map = ecm.get_dof_map();
  std::vector<unsigned int> dof_indices_ecm;

  // MDE system
  auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
  const unsigned int v_mde = mde.variable_number("mde");
  const DofMap &mde_map = mde.get_dof_map();
  std::vector<unsigned int> dof_indices_mde;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = mde.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfvfe::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_mde.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number mde_cur = 0.;
      Number mde_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Number ecm_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
        mde_old += phi[l][qp] * mde.old_solution(dof_indices_mde[l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
      }

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number nec_proj = util::project_concentration(nec_cur);
      Number ecm_proj = util::project_concentration(ecm_cur);
      Number mde_proj = util::project_concentration(mde_cur);

      Number aux_1 = deck->d_lambda_MDE_P * (tum_proj - nec_proj) *
                     ecm_proj * deck->d_sigma_HP / (1. + nut_proj);

      Number compute_rhs =
          JxW[qp] * (mde_old + dt * aux_1 * (1. - mde_proj) -
                     dt * deck->d_lambda_MDE_D * mde_proj -
                     dt * deck->d_lambda_ECM_D * ecm_proj * mde_proj);

      Number compute_mat = JxW[qp];

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_MDE * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    mde_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_mde);
    mde.matrix->add_matrix(Ke, dof_indices_mde);
    mde.rhs->add_vector(Fi, dof_indices_mde);
  }
}