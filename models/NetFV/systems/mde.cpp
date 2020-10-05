////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfv::initial_condition_mde(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "MDE");

  if (var_name == "mde") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    double val = 0.;
    for (unsigned int i = 0; i < deck->d_tum_ic_data.size(); i++)
      val += initial_condition_hyp_kernel(p, deck->d_dim,
                                          deck->d_tum_ic_data[i].d_ic_type,
                                          deck->d_tum_ic_data[i].d_ic_center,
                                          deck->d_tum_ic_data[i].d_tum_ic_radius,
                                          deck->d_tum_ic_data[i].d_hyp_ic_radius);

    return deck->d_mde_ic_val * val;
  }

  return 0.;
}

// Assembly class
void netfv::MdeAssembly::assemble() {
  assemble_face();
  assemble_1();
}

void netfv::MdeAssembly::assemble_face() {

  // call diffusion-advection calculation function
  if (d_model_p->get_input_deck().d_advection_active)
    netfv::assemble_diffusion_advection(*this,
                                        d_model_p->get_pres_assembly(),
                                        d_model_p->get_tum_assembly(), this->d_model_p);
  else
    netfv::assemble_diffusion(*this, this->d_model_p);
}

void netfv::MdeAssembly::assemble_1() {

  // Get required system alias
  // auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real mde_old = 0.;
  Real nut_cur = 0.;
  Real tum_cur = 0.;
  Real nec_cur = 0.;
  Real ecm_cur = 0.;

  Real nut_proj = 0.;
  Real tum_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // get fields at this element
    nut_cur = nut.get_current_sol(0);
    tum_cur = tum.get_current_sol_var(0, 0);
    nec_cur = nec.get_current_sol(0);
    ecm_cur = ecm.get_current_sol(0);
    mde_old = get_old_sol(0);

    if (deck.d_assembly_method == 1) {

      Real aux_1 = deck.d_lambda_MDE_P * (tum_cur - nec_cur) * ecm_cur *
                   deck.d_sigma_HP / (1. + nut_cur);

      compute_rhs = deck.d_elem_size * (mde_old + dt * aux_1);

      compute_mat =
        deck.d_elem_size * (1. + dt * deck.d_lambda_MDE_D + dt * aux_1 +
                            dt * deck.d_lambda_ECM_D * ecm_cur);
    } else {

      tum_proj = util::project_concentration(tum_cur);
      nec_proj = util::project_concentration(nec_cur);
      ecm_proj = util::project_concentration(ecm_cur);
      nut_proj = util::project_concentration(nut_cur);

      Real aux_1 = deck.d_lambda_MDE_P * (tum_proj - nec_proj) *
                   ecm_proj * deck.d_sigma_HP / (1. + nut_proj);

      compute_rhs = deck.d_elem_size * (mde_old + dt * aux_1);

      compute_mat = deck.d_elem_size * (1. + dt * deck.d_lambda_MDE_D + dt * aux_1 +
                                        dt * deck.d_lambda_ECM_D * ecm_proj);
    }

    // add
    Ke(0, 0) += compute_mat;

    // previous time step term
    Fe(0) += compute_rhs;

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}