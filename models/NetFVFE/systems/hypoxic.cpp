////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"

Number netfvfe::initial_condition_hyp_kernel(const Point &p,
    const netfvfe::InputDeck *deck) {

  const unsigned int dim = deck->d_dim;
  const unsigned int num_ic = deck->d_tum_ic_data.size();
  if (num_ic == 0)
    return 0.;

  for (unsigned int ic = 0; ic < num_ic; ic++) {

    auto data = deck->d_tum_ic_data[ic];

    const std::string type = data.d_ic_type;
    const Point xc =
        Point(data.d_ic_center[0], data.d_ic_center[1], data.d_ic_center[2]);
    const Point dx = p - xc;

    if (type == "tumor_hypoxic_spherical") {

      if (dx.norm() < data.d_tum_ic_radius[0] - 1.0E-12) {

        return 1. - util::exp_decay_function(
            dx.norm() / data.d_tum_ic_radius[0], 4.);

      } else if (dx.norm() > data.d_tum_ic_radius[0] - 1.0E-12 and
                 dx.norm() < data.d_hyp_ic_radius[0] - 1.0E-12) {

        return util::exp_decay_function(
            (dx.norm() - data.d_tum_ic_radius[0]) /
            (data.d_hyp_ic_radius[0] - data.d_tum_ic_radius[0]),
            4.);
      }
    } else if (type == "tumor_hypoxic_elliptical") {

      // transform ellipse into ball of radius
      double small_ball_r = 0.;
      for (unsigned int i = 0; i < dim; i++)
        small_ball_r = data.d_tum_ic_radius[i] * data.d_tum_ic_radius[i];
      small_ball_r = std::sqrt(small_ball_r);

      Point p_small_ball = util::ellipse_to_ball(p, xc, data.d_tum_ic_radius,
                                                 dim, small_ball_r);

      // transform ellipse into ball of radius
      double large_ball_r = 0.;
      for (unsigned int i = 0; i < dim; i++)
        large_ball_r = data.d_hyp_ic_radius[i] * data.d_hyp_ic_radius[i];
      large_ball_r = std::sqrt(large_ball_r);

      Point p_large_ball = util::ellipse_to_ball(p, xc, data.d_hyp_ic_radius,
                                                 dim, large_ball_r);

      if (p_small_ball.norm() < small_ball_r - 1.0E-12) {

        return 1. - util::exp_decay_function(
            p_small_ball.norm() / small_ball_r, 4.);

      } else if (p_small_ball.norm() > small_ball_r - 1.0E-12 and
                 p_small_ball.norm() < large_ball_r - 1.0E-12) {

        return util::exp_decay_function((p_small_ball.norm() - small_ball_r) /
                                        (large_ball_r - small_ball_r),
                                        4.);
      }
    }
  }

  return 0.;
}

Number netfvfe::initial_condition_hyp(const Point &p, const Parameters &es,
                                     const std::string &system_name,
                                     const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Hypoxic");

  if (var_name == "hypoxic") {

    const auto *deck = es.get<netfvfe::InputDeck *>("input_deck");

    return initial_condition_hyp_kernel(p, deck);
  }

  return 0.;
}

void netfvfe::HypAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfvfe::HypAssembly::assemble_1() {

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
  FEType fe_type = hyp.variable_type(0);
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
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_hyp.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number hyp_cur = 0.;
      Number hyp_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Gradient che_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        hyp_old += phi[l][qp] * hyp.old_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        che_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_tum_var[1][l]));
      }

      Number mobility =
          deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      Number compute_rhs =
          JxW[qp] * (hyp_old + dt * deck->d_lambda_PH *
                                   util::heaviside(deck->d_sigma_PH - nut_cur) *
                                   (tum_cur - nec_cur));

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_A +
                     dt * deck->d_lambda_HP *
                         util::heaviside(nut_cur - deck->d_sigma_HP) +
                     dt * deck->d_lambda_PH *
                         util::heaviside(deck->d_sigma_PH - nut_cur) +
                     dt * deck->d_lambda_HN *
                         util::heaviside(deck->d_sigma_HN - nut_cur));

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        // gradient term
        Fi(i) -= JxW[qp] * dt * mobility * che_grad * dphi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    hyp_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_hyp);
    hyp.matrix->add_matrix(Ke, dof_indices_hyp);
    hyp.rhs->add_vector(Fi, dof_indices_hyp);
  }
}

void netfvfe::HypAssembly::assemble_2() {

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
  FEType fe_type = hyp.variable_type(0);
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
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_hyp.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number hyp_cur = 0.;
      Number hyp_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Gradient che_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        hyp_old += phi[l][qp] * hyp.old_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        che_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_tum_var[1][l]));
      }

      Number mobility =
          deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number nec_proj = util::project_concentration(nec_cur);

      Number compute_rhs =
          JxW[qp] *
          (hyp_old + dt * deck->d_lambda_PH *
                         util::heaviside(deck->d_sigma_PH - nut_proj) *
                         (tum_proj - nec_proj));

      Number compute_mat =
          JxW[qp] * (1. + dt * deck->d_lambda_A +
                     dt * deck->d_lambda_HP *
                         util::heaviside(nut_proj - deck->d_sigma_HP) +
                     dt * deck->d_lambda_PH *
                         util::heaviside(deck->d_sigma_PH - nut_proj) +
                     dt * deck->d_lambda_HN *
                         util::heaviside(deck->d_sigma_HN - nut_proj));

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        // gradient term
        Fi(i) -= JxW[qp] * dt * mobility * che_grad * dphi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    hyp_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_hyp);
    hyp.matrix->add_matrix(Ke, dof_indices_hyp);
    hyp.rhs->add_vector(Fi, dof_indices_hyp);
  }
}

void netfvfe::HypAssembly::assemble_3() {

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
  FEType fe_type = hyp.variable_type(0);
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
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_hyp.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number hyp_cur = 0.;
      Number hyp_old = 0.;
      Number tum_cur = 0.;
      Number nec_cur = 0.;
      Gradient che_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        hyp_old += phi[l][qp] * hyp.old_solution(dof_indices_hyp[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        che_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_tum_var[1][l]));
      }

      Number mobility =
          deck->d_bar_M_H * pow(hyp_cur, 2) * pow(1. - hyp_cur, 2);

      // get projected values of species
      Number tum_proj = util::project_concentration(tum_cur);
      Number nec_proj = util::project_concentration(nec_cur);
      Number hyp_proj = util::project_concentration(hyp_cur);

      Number compute_rhs =
          JxW[qp] *
          (hyp_old +
           dt * deck->d_lambda_PH *
               util::heaviside(deck->d_sigma_PH - nut_proj) *
               (tum_proj - hyp_proj - nec_proj) -
           dt * deck->d_lambda_A * hyp_proj -
           dt * deck->d_lambda_HP *
               util::heaviside(nut_proj - deck->d_sigma_HP) * hyp_proj -
           dt * deck->d_lambda_HN *
               util::heaviside(deck->d_sigma_HN - nut_proj) * hyp_proj);

      Number compute_mat = JxW[qp];

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        // gradient term
        Fi(i) -= JxW[qp] * dt * mobility * che_grad * dphi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    hyp_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_hyp);
    hyp.matrix->add_matrix(Ke, dof_indices_hyp);
    hyp.rhs->add_vector(Fi, dof_indices_hyp);
  }
}