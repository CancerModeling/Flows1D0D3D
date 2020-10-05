////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"


Number netfc::initial_condition_tum(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Tumor");

  if (var_name == "tumor") {

    const auto *deck = es.get<netfc::InputDeck *>("input_deck");

    const unsigned int dim = deck->d_dim;
    const unsigned int num_ic = deck->d_tum_ic_data.size();
    if (num_ic == 0)
      return 0.;

    for (unsigned int ic = 0; ic < num_ic; ic++) {

      auto data = deck->d_tum_ic_data[ic];

      const std::string type = data.d_ic_type;
      const Point xc = Point(data.d_ic_center[0], data.d_ic_center[1],
                             data.d_ic_center[2]);
      const Point dx = p - xc;

      if (type == "tumor_spherical" or type == "tumor_hypoxic_spherical") {
        if (dx.norm() < data.d_tum_ic_radius[0] - 1.0E-12) {

          // out << "here tum ic\n";

          return util::exp_decay_function(dx.norm() / data.d_tum_ic_radius[0],
                                          4.);
        }
      } else if (type == "tumor_elliptical" or
                 type == "tumor_hypoxic_elliptical") {

        // transform ellipse into ball of radius
        double small_ball_r = 0.;
        for (unsigned int i = 0; i < dim; i++)
          small_ball_r = data.d_tum_ic_radius[i] * data.d_tum_ic_radius[i];
        small_ball_r = std::sqrt(small_ball_r);

        Point p_small_ball = util::ellipse_to_ball(p, xc, data.d_tum_ic_radius,
                                                   dim, small_ball_r);

        if (p_small_ball.norm() < small_ball_r - 1.0E-12) {

          return util::exp_decay_function(p_small_ball.norm() / small_ball_r,
                                          4.);
        }
      }
    }

    return 0.;
  }

  return 0.;
}

void netfc::TumAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfc::TumAssembly::assemble_1() {

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
  auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = tum.variable_type(0);
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
  DenseSubMatrix<Number> Ke_var[2][2] = {
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_tum.size();
    const unsigned int n_var_dofs = dof_indices_tum_var[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2; var_i++)
      for (unsigned int var_j = 0; var_j < 2; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number tum_cur = 0.;
      Number tum_old = 0.;
      Number hyp_cur = 0.;
      Number nec_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_old += phi[l][qp] * tum.old_solution(dof_indices_tum_var[0][l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
      }

      Number pro_cur = tum_cur - hyp_cur - nec_cur;
      Number mobility =
        deck->d_bar_M_P * pow(pro_cur, 2) * pow(1. - pro_cur, 2) +
        deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      // compute rhs
      Number compute_rhs_tum = JxW[qp] * (tum_old +
                                          dt * deck->d_lambda_P * nut_cur *
                                            pro_cur * (1. - tum_cur) +
                                          dt * deck->d_lambda_A * nec_cur);

      Number compute_rhs_mu = JxW[qp] * (deck->d_bar_E_phi_T * tum_old *
                                           (4.0 * pow(tum_old, 2) - 6.0 * tum_old - 1.) -
                                         deck->d_chi_c * nut_cur);

      // compute matrix
      Number compute_mat_tum = JxW[qp] * (1. + dt * deck->d_lambda_A);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        //-- Tumor --//
        Fe_var[0](i) += compute_rhs_tum * phi[i][qp];

        //-- Chemical Potential --//
        Fe_var[1](i) += compute_rhs_mu * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          //-- Tumor --//
          Ke_var[0][0](i, j) += compute_mat_tum * phi[j][qp] * phi[i][qp];

          // coupling with chemical potential
          Ke_var[0][1](i, j) +=
            JxW[qp] * dt * mobility * dphi[j][qp] * dphi[i][qp];

          //-- Chemical_tumor --//
          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          // coupling with tumor
          Ke_var[1][0](i, j) -= JxW[qp] * 3.0 * deck->d_bar_E_phi_T * phi[j][qp] * phi[i][qp];

          Ke_var[1][0](i, j) -=
            JxW[qp] * pow(deck->d_epsilon_T, 2) * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    tum_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                               dof_indices_tum);
    tum.matrix->add_matrix(Ke, dof_indices_tum);
    tum.rhs->add_vector(Fe, dof_indices_tum);
  }
}

void netfc::TumAssembly::assemble_2() {

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
  auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = tum.variable_type(0);
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
  DenseSubMatrix<Number> Ke_var[2][2] = {
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_tum.size();
    const unsigned int n_var_dofs = dof_indices_tum_var[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2; var_i++)
      for (unsigned int var_j = 0; var_j < 2; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number tum_cur = 0.;
      Number tum_old = 0.;
      Number hyp_cur = 0.;
      Number nec_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_old += phi[l][qp] * tum.old_solution(dof_indices_tum_var[0][l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
      }

      Number pro_cur = tum_cur - hyp_cur - nec_cur;
      Number mobility =
        deck->d_bar_M_P * pow(pro_cur, 2) * pow(1. - pro_cur, 2) +
        deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number tum_proj_old = util::project_concentration(tum_old);
      Number hyp_proj = util::project_concentration(hyp_cur);
      Number nec_proj = util::project_concentration(nec_cur);
      Number pro_proj = util::project_concentration(pro_cur);

      // compute rhs
      Number compute_rhs_tum = JxW[qp] * (tum_old +
                                          dt * deck->d_lambda_P * nut_proj *
                                            (tum_proj - hyp_proj - nec_proj) * (1. - tum_proj) +
                                          dt * deck->d_lambda_A * nec_proj);

      Number compute_rhs_mu =
        JxW[qp] *
        (deck->d_bar_E_phi_T * tum_proj_old *
           (4.0 * pow(tum_proj_old, 2) - 6.0 * tum_proj_old - 1.) -
         deck->d_chi_c * nut_proj);

      // compute matrix
      Number compute_mat_tum = JxW[qp] * (1. + dt * deck->d_lambda_A);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        //-- Tumor --//
        Fe_var[0](i) += compute_rhs_tum * phi[i][qp];

        //-- Chemical Potential --//
        Fe_var[1](i) += compute_rhs_mu * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          //-- Tumor --//
          Ke_var[0][0](i, j) += compute_mat_tum * phi[j][qp] * phi[i][qp];

          // coupling with chemical potential
          Ke_var[0][1](i, j) +=
            JxW[qp] * dt * mobility * dphi[j][qp] * dphi[i][qp];

          //-- Chemical_tumor --//
          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          // coupling with tumor
          Ke_var[1][0](i, j) -= JxW[qp] * 3.0 * deck->d_bar_E_phi_T * phi[j][qp] * phi[i][qp];

          Ke_var[1][0](i, j) -=
            JxW[qp] * pow(deck->d_epsilon_T, 2) * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    tum_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                               dof_indices_tum);
    tum.matrix->add_matrix(Ke, dof_indices_tum);
    tum.rhs->add_vector(Fe, dof_indices_tum);
  }
}

void netfc::TumAssembly::assemble_3() {

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
  auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = tum.variable_type(0);
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
  DenseSubMatrix<Number> Ke_var[2][2] = {
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
    {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_tum.size();
    const unsigned int n_var_dofs = dof_indices_tum_var[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2; var_i++)
      for (unsigned int var_j = 0; var_j < 2; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number tum_cur = 0.;
      Number tum_old = 0.;
      Number hyp_cur = 0.;
      Number nec_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_old += phi[l][qp] * tum.old_solution(dof_indices_tum_var[0][l]);
        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
      }

      Number pro_cur = tum_cur - hyp_cur - nec_cur;
      Number mobility =
        deck->d_bar_M_P * pow(pro_cur, 2) * pow(1. - pro_cur, 2) +
        deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number tum_proj_old = util::project_concentration(tum_old);
      Number hyp_proj = util::project_concentration(hyp_cur);
      Number nec_proj = util::project_concentration(nec_cur);
      Number pro_proj = util::project_concentration(pro_cur);

      // compute rhs
      Number compute_rhs_tum = JxW[qp] * (tum_old +
                                          dt * deck->d_lambda_P * nut_proj *
                                            (tum_proj - hyp_proj - nec_proj) *
                                            (1. - tum_proj) -
                                          dt * deck->d_lambda_A *
                                            (tum_proj - nec_proj));

      Number compute_rhs_mu =
        JxW[qp] *
        (deck->d_bar_E_phi_T * tum_proj_old *
           (4.0 * pow(tum_proj_old, 2) - 6.0 * tum_proj_old - 1.) -
         deck->d_chi_c * nut_proj);

      // compute matrix
      Number compute_mat_tum = JxW[qp];

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        //-- Tumor --//
        Fe_var[0](i) += compute_rhs_tum * phi[i][qp];

        //-- Chemical Potential --//
        Fe_var[1](i) += compute_rhs_mu * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          //-- Tumor --//
          Ke_var[0][0](i, j) += compute_mat_tum * phi[j][qp] * phi[i][qp];

          // coupling with chemical potential
          Ke_var[0][1](i, j) +=
            JxW[qp] * dt * mobility * dphi[j][qp] * dphi[i][qp];

          //-- Chemical_tumor --//
          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          // coupling with tumor
          Ke_var[1][0](i, j) -= JxW[qp] * 3.0 * deck->d_bar_E_phi_T * phi[j][qp] * phi[i][qp];

          Ke_var[1][0](i, j) -=
            JxW[qp] * pow(deck->d_epsilon_T, 2) * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    tum_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                               dof_indices_tum);
    tum.matrix->add_matrix(Ke, dof_indices_tum);
    tum.rhs->add_vector(Fe, dof_indices_tum);
  }
}