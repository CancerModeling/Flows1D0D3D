////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "systems.hpp"
#include "../model.hpp"

Number netfc::initial_condition_ecm(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"ECM");

  if (var_name == "ecm") {

    const auto *deck = es.get<netfc::InputDeck *>("input_deck");
    const auto &ic_data = deck->d_ecm_ic_data;

    if (ic_data.d_type.empty())
      return 0.;
    else if (ic_data.d_type == "spherical") {

      Point dx = p - Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
          ic_data.d_geom_params[2]);
      double r = ic_data.d_geom_params[3];

      if (dx.norm() < r - 1.0E-12)
        return ic_data.d_val * util::exp_decay_function(dx.norm() / r, 4.);
      else
        return 0.;
    } else if (ic_data.d_type == "elliptical") {

      Point xc = Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
                       ic_data.d_geom_params[2]);
      std::vector<double> r = {ic_data.d_geom_params[3], ic_data
                               .d_geom_params[4], ic_data.d_geom_params[5]};
      const Point dx = p - xc;

      // transform ellipse into ball of radius
      double ball_r = 0.;
      for (unsigned int i = 0; i < deck->d_dim; i++)
        ball_r = r[i] * r[i];
      ball_r = std::sqrt(ball_r);

      Point p_ball = util::ellipse_to_ball(p, xc, r,
                                                 deck->d_dim, ball_r);

      if (p_ball.norm() < ball_r - 1.0E-12) {

        return ic_data.d_val * util::exp_decay_function(p_ball.norm() / ball_r, 4.);
      } else
        return 0.;

    } else if (ic_data.d_type == "box") {

      Point x1 = Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
                           ic_data.d_geom_params[2]);
      Point x2 = Point(ic_data.d_geom_params[3],ic_data.d_geom_params[4],
                       ic_data.d_geom_params[5]);

      if (util::is_inside_box(p, {x1, x2}))
        return ic_data.d_val;
      else
        return 0.;
    } else if (ic_data.d_type == "constant") {
      return ic_data.d_val;
    }
  }

  return 0.;
}

void netfc::boundary_condition_ecm(EquationSystems &es) {

  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");
}

void netfc::EcmAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfc::EcmAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

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
  FEType fe_type = ecm.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_ecm.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number ecm_cur = 0.;
      Number ecm_old = 0.;
      Number mde_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
        ecm_old += phi[l][qp] * ecm.old_solution(dof_indices_ecm[l]);
        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
      }

      Number compute_rhs = JxW[qp] * (ecm_old + dt * deck->d_lambda_ECM_P *
          nut_cur * util::heaviside(ecm_cur - deck->d_bar_phi_ECM_P));

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_ECM_D * mde_cur +
                     dt * deck->d_lambda_ECM_P * nut_cur *
                         util::heaviside(ecm_cur - deck->d_bar_phi_ECM_P));

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    ecm_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_ecm);
    ecm.matrix->add_matrix(Ke, dof_indices_ecm);
    ecm.rhs->add_vector(Fi, dof_indices_ecm);
  }
}

void netfc::EcmAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

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
  FEType fe_type = ecm.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_ecm.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number ecm_cur = 0.;
      Number ecm_old = 0.;
      Number mde_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
        ecm_old += phi[l][qp] * ecm.old_solution(dof_indices_ecm[l]);
        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
      }

      // get projected values of species
      Number ecm_proj = util::project_concentration(ecm_cur);
      Number mde_proj = util::project_concentration(mde_cur);

      Number compute_rhs = JxW[qp] * (ecm_old + dt * deck->d_lambda_ECM_P *
                                                nut_proj * util::heaviside(ecm_proj - deck->d_bar_phi_ECM_P));

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_ECM_D * mde_proj +
                     dt * deck->d_lambda_ECM_P * nut_proj *
                     util::heaviside(ecm_proj - deck->d_bar_phi_ECM_P));

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    ecm_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_ecm);
    ecm.matrix->add_matrix(Ke, dof_indices_ecm);
    ecm.rhs->add_vector(Fi, dof_indices_ecm);
  }
}

void netfc::EcmAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Nutrient system
  auto &nut =
      es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

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
  FEType fe_type = ecm.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Model parameters
  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    const unsigned int n_dofs = dof_indices_ecm.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number ecm_cur = 0.;
      Number ecm_old = 0.;
      Number mde_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
        ecm_old += phi[l][qp] * ecm.old_solution(dof_indices_ecm[l]);
        mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
      }

      // get projected values of species
      Number ecm_proj = util::project_concentration(ecm_cur);
      Number mde_proj = util::project_concentration(mde_cur);

      Number compute_rhs =
          JxW[qp] *
          (ecm_old + dt * deck->d_lambda_ECM_P * nut_proj * (1. - ecm_proj) *
                         util::heaviside(ecm_proj - deck->d_bar_phi_ECM_P)
           - dt * deck->d_lambda_ECM_D * mde_proj * ecm_proj);

      Number compute_mat = JxW[qp];

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    ecm_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_ecm);
    ecm.matrix->add_matrix(Ke, dof_indices_ecm);
    ecm.rhs->add_vector(Fi, dof_indices_ecm);
  }
}