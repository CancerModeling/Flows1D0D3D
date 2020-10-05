////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "systems.hpp"
#include "../model.hpp"

Number netfc::initial_condition_nec(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Necrotic");

  return 0.;
}

void netfc::NecAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2 || deck.d_assembly_method == 3)
    assemble_2();
}

void netfc::NecAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum =
    es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(2);

  // Nutrient system
  auto &nut =
    es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");

  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp =
    es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec =
    es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = nec.variable_type(0);
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

    tum_map.dof_indices(elem, dof_indices);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_nec.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number nec_old = 0.;
      Number nec_cur = 0.;
      Number hyp_cur = 0.;
      Gradient che_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        nec_old += phi[l][qp] * nec.old_solution(dof_indices_nec[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        che_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_var[1][l]));
      }

      Number compute_rhs = JxW[qp] *
                           (nec_old + dt * deck->d_lambda_HN *
                                        util::heaviside(deck->d_sigma_HN - nut_cur) *
                                        hyp_cur);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    nec_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_nec);
    nec.matrix->add_matrix(Ke, dof_indices_nec);
    nec.rhs->add_vector(Fi, dof_indices_nec);
  }
}

void netfc::NecAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum =
    es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(2);

  // Nutrient system
  auto &nut =
    es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");

  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp =
    es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec =
    es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = nec.variable_type(0);
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

    tum_map.dof_indices(elem, dof_indices);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);

    const unsigned int n_dofs = dof_indices_nec.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    // get nutrient at this element
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    Number nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number nec_old = 0.;
      Number nec_cur = 0.;
      Number hyp_cur = 0.;
      Gradient che_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        nec_old += phi[l][qp] * nec.old_solution(dof_indices_nec[l]);
        nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
        hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
        che_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_var[1][l]));
      }

      // get projected values of species
      Number hyp_proj = util::project_concentration(hyp_cur);

      Number compute_rhs =
        JxW[qp] *
        (nec_old + dt * deck->d_lambda_HN *
                     util::heaviside(deck->d_sigma_HN - nut_proj) *
                     hyp_proj);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += compute_rhs * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    nec_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fi,
                                                               dof_indices_nec);
    nec.matrix->add_matrix(Ke, dof_indices_nec);
    nec.rhs->add_vector(Fi, dof_indices_nec);
  }
}
