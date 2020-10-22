////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfvfe::initial_condition_hyp(const Point &p, const Parameters &es,
                                      const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Hypoxic");

  if (var_name == "hypoxic") {
    return 0.;
  }
}

// Assembly class
void netfvfe::HypAssembly::assemble() {

  assemble_1();
}

void netfvfe::HypAssembly::assemble_1() {

  // Get required system alias
  auto &nut = d_model_p->get_nut_assembly();
  auto &pro = d_model_p->get_pro_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();
  auto &vel = d_model_p->get_vel_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real advection_factor = deck.d_advection_active ? 1. : 0.;

  // Store current and old solution
  Real tum_cur = 0.;
  Real nut_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;
  Real pro_cur = 0.;
  Real ecm_cur = 0.;
  Real pro_old = 0.;
  Real tum_old = 0.;
  Real hyp_old = 0.;
  Real nec_old = 0.;
  Real ecm_old = 0.;

  Real mobility = 0.;

  Gradient vel_cur = 0.;

  Real compute_rhs_hyp = 0.;
  Real compute_rhs_mu = 0.;
  Real compute_mat_hyp = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    hyp.init_dof(elem);
    nut.init_dof(elem);
    pro.init_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);
    vel.init_dof(elem);

    // init fe and element matrix and vector
    hyp.init_fe(elem);

    // get finite-volume quantities
    nut_cur = nut.get_current_sol(0);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      if (d_implicit_assembly) {
        // require
        // old: pro, hyp, nec
        // new: nut, vel, pro, hyp, nec, ecm
        pro_old = 0.;
        hyp_old = 0.;
        nec_old = 0.;
        tum_old = 0.;
        vel_cur = 0.;
        pro_cur = 0.;
        hyp_cur = 0.;
        nec_cur = 0.;
        ecm_cur = 0.;

        for (unsigned int l = 0; l < d_phi.size(); l++) {

          pro_old += d_phi[l][qp] * pro.get_old_sol_var(l, 0);
          hyp_old += d_phi[l][qp] * hyp.get_old_sol_var(l, 0);
          nec_old += d_phi[l][qp] * nec.get_old_sol(l);

          pro_cur += d_phi[l][qp] * pro.get_current_sol_var(l, 0);
          hyp_cur += d_phi[l][qp] * hyp.get_current_sol_var(l, 0);
          nec_cur += d_phi[l][qp] * nec.get_current_sol(l);
          ecm_cur += d_phi[l][qp] * ecm.get_current_sol(l);

          for (unsigned int ll = 0; ll < d_mesh.mesh_dimension(); ll++)
            vel_cur(ll) += d_phi[l][qp] * vel.get_current_sol_var(l, ll);
        }

        tum_old = pro_old + hyp_old + nec_old;
        tum_cur = pro_cur + hyp_cur + nec_cur;

        mobility = deck.d_bar_M_H * pow(util::proj(hyp_cur) * util::proj(1. - tum_cur), 2);

        compute_rhs_hyp =
          d_JxW[qp] * (hyp_old + dt * deck.d_lambda_PH * util::heaviside(deck.d_sigma_PH - nut_cur) * pro_cur);

        compute_rhs_mu =
          d_JxW[qp] * (deck.d_bar_E_phi_T * tum_old *
                         (4.0 * pow(tum_old, 2) - 6.0 * tum_old - 1.) +
                       deck.d_bar_E_phi_H * hyp_old *
                         (4.0 * pow(hyp_old, 2) - 6.0 * hyp_old - 1.) +
                       3. * deck.d_bar_E_phi_T * (pro_cur + nec_cur) -
                       deck.d_chi_c * nut_cur - deck.d_chi_h * ecm_cur);

        compute_mat_hyp =
          d_JxW[qp] * (1. + dt * deck.d_lambda_A -
                       dt * deck.d_lambda_Ph * nut_cur * util::proj(1. - tum_cur) +
                       dt * deck.d_lambda_HP *
                         util::heaviside(nut_cur - deck.d_sigma_HP) +
                       dt * deck.d_lambda_HN *
                         util::heaviside(deck.d_sigma_HN - nut_cur));

      } else {
        // require
        // old: pro, hyp, nec, ecm
        // new: nut, vel, pro
        pro_old = 0.;
        hyp_old = 0.;
        nec_old = 0.;
        ecm_old = 0.;
        tum_old = 0.;
        vel_cur = 0.;
        pro_cur = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {

          pro_old += d_phi[l][qp] * pro.get_old_sol_var(l, 0);
          hyp_old += d_phi[l][qp] * hyp.get_old_sol_var(l, 0);
          nec_old += d_phi[l][qp] * nec.get_old_sol(l);
          ecm_old += d_phi[l][qp] * ecm.get_old_sol(l);

          pro_cur += d_phi[l][qp] * pro.get_current_sol_var(l, 0);

          for (unsigned int ll = 0; ll < d_mesh.mesh_dimension(); ll++)
            vel_cur(ll) += d_phi[l][qp] * vel.get_current_sol_var(l, ll);
        }

        tum_old = pro_old + hyp_old + nec_old;

        mobility = deck.d_bar_M_P * pow(util::proj(hyp_old) * util::proj(1. - tum_old), 2);

        compute_rhs_hyp =
          d_JxW[qp] * (hyp_old + dt * deck.d_lambda_PH * util::heaviside(deck.d_sigma_PH - nut_cur) * pro_cur - dt * deck.d_lambda_HP * util::heaviside(nut_cur - deck.d_sigma_HP) * hyp_old);

        compute_rhs_mu =
          d_JxW[qp] * (deck.d_bar_E_phi_T * tum_old *
                         (4.0 * pow(tum_old, 2) - 6.0 * tum_old - 1.) +
                       deck.d_bar_E_phi_H * hyp_old *
                         (4.0 * pow(hyp_old, 2) - 6.0 * hyp_old - 1.) +
                       3. * deck.d_bar_E_phi_T * (pro_old + nec_old) -
                       deck.d_chi_c * nut_cur - deck.d_chi_h * ecm_old);

        compute_mat_hyp =
          d_JxW[qp] * (1. + dt * deck.d_lambda_A -
                       dt * deck.d_lambda_Ph * nut_cur * util::proj(1. - tum_old) +
                       dt * deck.d_lambda_HN *
                         util::heaviside(deck.d_sigma_HN - nut_cur));
      }

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        // hypoxic
        d_Fe_var[0](i) += compute_rhs_hyp * d_phi[i][qp];

        // chemical
        d_Fe_var[1](i) += compute_rhs_mu * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          // hypoxic
          d_Ke_var[0][0](i, j) += compute_mat_hyp * d_phi[j][qp] *
                                  d_phi[i][qp];

          // advection of hypoxic
          d_Ke_var[0][0](i, j) -=
            advection_factor * d_JxW[qp] * dt * d_phi[j][qp] * vel_cur * d_dphi[i][qp];

          // coupling with chemical potential
          d_Ke_var[0][1](i, j) +=
            d_JxW[qp] * dt * mobility * d_dphi[j][qp] * d_dphi[i][qp];

          // chemical
          d_Ke_var[1][1](i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];

          // coupling with tumor
          d_Ke_var[1][0](i, j) -= d_JxW[qp] * 3.0 * (deck.d_bar_E_phi_T + deck.d_bar_E_phi_H) * d_phi[j][qp] * d_phi[i][qp];

          d_Ke_var[1][0](i, j) -= d_JxW[qp] * pow(deck.d_epsilon_H, 2) *
                                  d_dphi[j][qp] * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}
