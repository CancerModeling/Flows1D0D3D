////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "systems.hpp"
#include "../model.hpp"

Number netfc::initial_condition_taf(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"TAF");

  return 0.;
}

Number netfc::initial_condition_grad_taf(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"TAF_Gradient");

  return 0.;
}

void netfc::TafAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfc::TafAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // TAF system
  auto &taf =
      es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // Hypoxic system
  auto &hyp =
      es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = taf.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
  const std::vector<Point> &q_points = fe->get_xyz();

  // Model parameters
  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    taf_map.dof_indices(elem, dof_indices_taf, v_taf);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);

    const unsigned int n_dofs = dof_indices_taf.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number taf_cur = 0.;
      Number taf_old = 0.;
      Number hyp_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        taf_cur += phi[l][qp] * taf.current_solution(dof_indices_taf[l]);
        taf_old += phi[l][qp] * taf.old_solution(dof_indices_taf[l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
      }

      Number compute_rhs = JxW[qp] * (taf_old + dt * deck->d_lambda_TAF *
          hyp_cur);

      // ADD ARTIFICIAL SOURCE
      {
        const Point x = q_points[qp];
        double L_half = 0.5 * deck->d_domain_params[1];
        const Point xc = Point(L_half, L_half, L_half);
        // spherical source
        if (false) {
          if ((x - xc).norm() < 0.1 * L_half)
            compute_rhs += JxW[qp] * dt * deck->d_lambda_TAF;
        }

        // source on cylinder along z axis
        {
          if (x(0) > 0.9 * L_half and x(0) < 1.1 * L_half and
              x(1) > 0.9 * L_half and x(1) < 1.1 * L_half) {

            //            compute_rhs += JxW[qp] * dt * deck->d_lambda_TAF *
            //                           std::sin(2. * (2. * M_PI) * x(2) / L_half);
            const Point x_plane = Point(x(0), x(1), 0.);
            const Point xc_plane = Point(L_half, L_half, 0.);
            if ((x_plane - xc_plane).norm() < 0.1 * L_half)
              compute_rhs +=
                  JxW[qp] * dt * deck->d_lambda_TAF;
//                  util::exp_decay_function((x_plane - xc_plane).norm() / (0.1 *
//                  L_half), 4);
          }
        }
      }

      Number compute_mat = JxW[qp] * (1. + dt * deck->d_lambda_TAF * hyp_cur);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_TAF * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    taf_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_taf);
    taf.matrix->add_matrix(Ke, dof_indices_taf);
    taf.rhs->add_vector(Fi, dof_indices_taf);
  }
}

void netfc::TafAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // TAF system
  auto &taf =
      es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // Hypoxic system
  auto &hyp =
      es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = taf.variable_type(0);
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

    taf_map.dof_indices(elem, dof_indices_taf, v_taf);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);

    const unsigned int n_dofs = dof_indices_taf.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number taf_cur = 0.;
      Number taf_old = 0.;
      Number hyp_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        taf_cur += phi[l][qp] * taf.current_solution(dof_indices_taf[l]);
        taf_old += phi[l][qp] * taf.old_solution(dof_indices_taf[l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
      }

      // get projected values of species
      Number hyp_proj = util::project_concentration(hyp_cur);
      Number taf_proj = util::project_concentration(taf_cur);

      Number compute_rhs = JxW[qp] * (taf_old + dt * deck->d_lambda_TAF *
                                                hyp_proj);

      Number compute_mat = JxW[qp] * (1. + dt * deck->d_lambda_TAF * hyp_proj);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_TAF * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    taf_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_taf);
    taf.matrix->add_matrix(Ke, dof_indices_taf);
    taf.rhs->add_vector(Fi, dof_indices_taf);
  }
}

void netfc::TafAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // TAF system
  auto &taf =
      es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // Hypoxic system
  auto &hyp =
      es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = taf.variable_type(0);
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

    taf_map.dof_indices(elem, dof_indices_taf, v_taf);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);

    const unsigned int n_dofs = dof_indices_taf.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number taf_cur = 0.;
      Number taf_old = 0.;
      Number hyp_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++) {

        taf_cur += phi[l][qp] * taf.current_solution(dof_indices_taf[l]);
        taf_old += phi[l][qp] * taf.old_solution(dof_indices_taf[l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
      }

      // get projected values of species
      Number hyp_proj = util::project_concentration(hyp_cur);
      Number taf_proj = util::project_concentration(taf_cur);

      Number compute_rhs =
          JxW[qp] *
          (taf_old + dt * deck->d_lambda_TAF * hyp_proj * (1. - taf_proj));

      Number compute_mat = JxW[qp];

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += compute_mat * phi[j][qp] * phi[i][qp];

          // gradient term
          Ke(i, j) +=
              JxW[qp] * dt * deck->d_D_TAF * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    taf_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_taf);
    taf.matrix->add_matrix(Ke, dof_indices_taf);
    taf.rhs->add_vector(Fi, dof_indices_taf);
  }
}

void netfc::GradTafAssembly::assemble() {
  const auto &deck = d_model_p->get_input_deck();
  if (deck.d_dim > 2)
    assemble_3d();
  else
    assemble_2d();
}

void netfc::GradTafAssembly::assemble_3d() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // TAF system
  auto &taf =
      es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // Gradient of TAF
  auto &grad =
      es.get_system<TransientLinearImplicitSystem>("TAF_Gradient");
  std::vector<unsigned int> v_grad(3);
  v_grad[0] = grad.variable_number("taf_gradx");
  v_grad[1] = grad.variable_number("taf_grady");
  v_grad[2] = grad.variable_number("taf_gradz");

  const DofMap &grad_map = grad.get_dof_map();
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(3);

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = grad.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[3][3] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[3] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    grad_map.dof_indices(elem, dof_indices);
    for (unsigned int var = 0; var < 3; var++)
      grad_map.dof_indices(elem, dof_indices_var[var], v_grad[var]);

    taf_map.dof_indices(elem, dof_indices_taf, v_taf);

    const unsigned int n_dofs = dof_indices.size();
    const unsigned int n_var_dofs = dof_indices_var[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 3; var_i++)
      for (unsigned int var_j = 0; var_j < 3; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 3; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Gradient taf_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        taf_grad.add_scaled(dphi[l][qp],
                            taf.current_solution(dof_indices_taf[l]));
      }

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fe_var[0](i) += JxW[qp] * taf_grad(0) * phi[i][qp];
        Fe_var[1](i) += JxW[qp] * taf_grad(1) * phi[i][qp];
        Fe_var[2](i) += JxW[qp] * taf_grad(2) * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke_var[0][0](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          Ke_var[2][2](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    grad_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                                dof_indices);
    grad.matrix->add_matrix(Ke, dof_indices);
    grad.rhs->add_vector(Fe, dof_indices);
  }
}

void netfc::GradTafAssembly::assemble_2d() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // TAF system
  auto &taf =
      es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // Gradient of TAF
  auto &grad =
      es.get_system<TransientLinearImplicitSystem>("TAF_Gradient");
  std::vector<unsigned int> v_grad(2);
  v_grad[0] = grad.variable_number("taf_gradx");
  v_grad[1] = grad.variable_number("taf_grady");

  const DofMap &grad_map = grad.get_dof_map();
  std::vector<unsigned int> dof_indices_grad;
  std::vector<std::vector<dof_id_type>> dof_indices_grad_var(2);

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = grad.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[2][2] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    grad_map.dof_indices(elem, dof_indices_grad);
    for (unsigned int var = 0; var < 2; var++)
      grad_map.dof_indices(elem, dof_indices_grad_var[var], v_grad[var]);

    taf_map.dof_indices(elem, dof_indices_taf, v_taf);

    const unsigned int n_dofs = dof_indices_grad.size();
    const unsigned int n_var_dofs = dof_indices_grad_var[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2; var_i++)
      for (unsigned int var_j = 0; var_j < 2; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Gradient taf_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        taf_grad.add_scaled(dphi[l][qp],
                            taf.current_solution(dof_indices_taf[l]));
      }

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fe_var[0](i) += JxW[qp] * taf_grad(0) * phi[i][qp];
        Fe_var[1](i) += JxW[qp] * taf_grad(1) * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke_var[0][0](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];

          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    grad_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                                dof_indices_grad);
    grad.matrix->add_matrix(Ke, dof_indices_grad);
    grad.rhs->add_vector(Fe, dof_indices_grad);
  }
}