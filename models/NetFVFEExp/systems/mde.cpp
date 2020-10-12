////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfvfeexp::initial_condition_mde(const Point &p, const Parameters &es,
                                         const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "MDE");

  if (var_name == "mde") {
    return 0.;
  }
}

// Assembly class
void netfvfeexp::MdeAssembly::assemble() {
  assemble_1();
}

void netfvfeexp::MdeAssembly::assemble_1() {

  // Get required system alias
  auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &pro = d_model_p->get_pro_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();
  auto &vel = d_model_p->get_vel_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real advection_factor = deck.d_advection_active ? 1. : 0.;

  // Store current and old solution
  Real mde_old = 0.;
  Real nut_cur = 0.;
  Real pro_cur = 0.;
  Real hyp_cur = 0.;
  Real ecm_cur = 0.;

  Real pro_old = 0.;
  Real hyp_old = 0.;
  Real ecm_old = 0.;

  Gradient vel_cur = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    mde.init_dof(elem);
    nut.init_dof(elem);
    pro.init_dof(elem);
    hyp.init_dof(elem);
    ecm.init_dof(elem);
    vel.init_dof(elem);

    // init fe and element matrix and vector
    mde.init_fe(elem);

    // get finite-volume quantities
    nut_cur = nut.get_current_sol(0);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      if (d_implicit_assembly) {
        // required
        // old: mde
        // new: nut, vel, pro, hyp, ecm
        mde_old = 0.;
        pro_cur = 0.;
        hyp_cur = 0.;
        ecm_cur = 0.;
        vel_cur = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {

          mde_old += d_phi[l][qp] * mde.get_old_sol(l);
          pro_cur += d_phi[l][qp] * pro.get_current_sol_var(l, 0);
          hyp_cur += d_phi[l][qp] * hyp.get_current_sol_var(l, 0);
          ecm_cur += d_phi[l][qp] * ecm.get_current_sol(l);

          for (unsigned int ll = 0; ll < d_mesh.mesh_dimension(); ll++)
            vel_cur(ll) += d_phi[l][qp] * vel.get_current_sol_var(l, ll);
        }

        Real aux_1 = deck.d_lambda_MDE_P * (pro_cur + hyp_cur) *
                     ecm_cur * deck.d_sigma_HP / (deck.d_sigma_HP + nut_cur);

        compute_rhs = d_JxW[qp] * (mde_old + dt * aux_1);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_MDE_D + dt * aux_1 +
                                   dt * deck.d_lambda_ECM_D * ecm_cur);

      } else {
        // required
        // old: mde, pro, hyp, ecm
        // new: nut, vel
        mde_old = 0.;
        pro_old = 0.;
        hyp_old = 0.;
        ecm_old = 0.;
        vel_cur = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {

          mde_old += d_phi[l][qp] * mde.get_old_sol(l);
          pro_old += d_phi[l][qp] * pro.get_old_sol_var(l, 0);
          hyp_old += d_phi[l][qp] * hyp.get_old_sol_var(l, 0);
          ecm_old += d_phi[l][qp] * ecm.get_old_sol(l);

          for (unsigned int ll = 0; ll < d_mesh.mesh_dimension(); ll++)
            vel_cur(ll) += d_phi[l][qp] * vel.get_current_sol_var(l, ll);
        }

        Real aux_1 = deck.d_lambda_MDE_P * (pro_old + hyp_old) *
                     ecm_old * deck.d_sigma_HP / (deck.d_sigma_HP + nut_cur);

        compute_rhs = d_JxW[qp] * (mde_old + dt * aux_1 - dt * deck.d_lambda_ECM_D * ecm_old);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_MDE_D + dt * aux_1);
      }

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe(i) += compute_rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke(i, j) += compute_mat * d_phi[j][qp] * d_phi[i][qp];

          // gradient term
          d_Ke(i, j) +=
            d_JxW[qp] * dt * deck.d_D_MDE * d_dphi[j][qp] * d_dphi[i][qp];

          // advection of mde
          d_Ke(i, j) -=
            advection_factor * d_JxW[qp] * dt * d_phi[j][qp] * vel_cur * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}
