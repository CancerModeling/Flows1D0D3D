////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfvfe::initial_condition_hyp(const Point &p, const Parameters &es,
                                      const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"Hypoxic");

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
  Real pro_old = 0.;
  Real tum_old = 0.;
  Real hyp_old = 0.;
  Real nec_old = 0.;
  Real ecm_cur = 0.;

  Real nut_proj = 0.;
  Real tum_proj = 0.;
  Real nec_proj = 0.;
  Real hyp_proj = 0.;
  Real pro_proj = 0.;
  Real ecm_proj = 0.;

  Real mobility = 0.;

  Gradient vel_cur = 0.;

  Real compute_rhs_hyp = 0.;
  Real compute_rhs_mu = 0.;
  Real compute_mat_hyp = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    pro.init_dof(elem);
    nec.init_dof(elem);
    ecm.init_dof(elem);
    vel.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    // get finite-volume quantities
    nut_cur = nut.get_current_sol(0);
    nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      pro_old = 0.; pro_cur = 0.; tum_cur = 0.; hyp_cur = 0.;
      nec_cur = 0.; ecm_cur = 0.;
      tum_old = 0.; hyp_old = 0.; nec_old = 0.;
      vel_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        hyp_old += d_phi[l][qp] * get_old_sol_var(l, 0);
        hyp_cur += d_phi[l][qp] * get_current_sol_var(l, 0);
        pro_cur += d_phi[l][qp] * pro.get_current_sol_var(l, 0);
        pro_old += d_phi[l][qp] * pro.get_old_sol_var(l, 0);
        nec_cur += d_phi[l][qp] * nec.get_current_sol(l);
        nec_old += d_phi[l][qp] * nec.get_old_sol(l);
        ecm_cur += d_phi[l][qp] * ecm.get_current_sol(l);

        for (unsigned int ll=0; ll<d_mesh.mesh_dimension(); ll++)
          vel_cur(ll) += d_phi[l][qp] * vel.get_current_sol_var(l, ll);
      }

      tum_cur = pro_cur + hyp_cur + nec_cur;
      tum_old = pro_old + hyp_old + nec_old;

      // get projected solution
      hyp_proj = util::project_concentration(hyp_cur);
      pro_proj = util::project_concentration(pro_cur);
      nec_proj = util::project_concentration(nec_cur);
      ecm_proj = util::project_concentration(ecm_cur);
      tum_proj = util::project_concentration(pro_proj + hyp_proj + nec_proj);

      mobility = deck.d_bar_M_H * pow(hyp_proj, 2) * pow(1. - hyp_proj, 2);

      if (deck.d_assembly_method == 1) {

        // compute quantities independent of dof loop
        compute_rhs_hyp =
            d_JxW[qp] * (hyp_old + dt * deck.d_lambda_PH * util::heaviside
                (deck.d_sigma_PH - nut_cur) * pro_cur);

        // keep the factor d_bar_E_phi_T in double well same for
        // prolific and hypoxic
        compute_rhs_mu =
            d_JxW[qp] * (deck.d_bar_E_phi_T * tum_old *
                         (4.0 * pow(tum_old, 2) - 6.0 * tum_old - 1.) +
                         3. * deck.d_bar_E_phi_T * (pro_cur + nec_cur) -
                         deck.d_chi_c * nut_cur - deck.d_chi_h * ecm_cur);

        compute_mat_hyp =
            d_JxW[qp] * (1. + dt * deck.d_lambda_A -
                         dt * deck.d_lambda_Ph * nut_cur +
                         dt * deck.d_lambda_HP *
                         util::heaviside(nut_cur - deck.d_sigma_HP) +
                         dt * deck.d_lambda_HN *
                         util::heaviside(deck.d_sigma_HN - nut_cur));
      } else {

        // compute quantities independent of dof loop
        compute_rhs_hyp =
            d_JxW[qp] * (hyp_old + dt * deck.d_lambda_PH * util::heaviside
                (deck.d_sigma_PH - nut_proj) * pro_proj);

        // keep the factor d_bar_E_phi_T in double well same for
        // prolific and hypoxic
        compute_rhs_mu =
            d_JxW[qp] * (deck.d_bar_E_phi_T * tum_old *
                         (4.0 * pow(tum_old, 2) - 6.0 * tum_old - 1.) +
                         3. * deck.d_bar_E_phi_T * (pro_proj + nec_proj) -
                         deck.d_chi_c * nut_proj - deck.d_chi_h * ecm_proj);

        compute_mat_hyp =
            d_JxW[qp] * (1. + dt * deck.d_lambda_A -
                         dt * deck.d_lambda_Ph * nut_proj +
                         dt * deck.d_lambda_HP *
                         util::heaviside(nut_proj - deck.d_sigma_HP) +
                         dt * deck.d_lambda_HN *
                         util::heaviside(deck.d_sigma_HN - nut_proj));
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
          d_Ke_var[1][0](i, j) -= d_JxW[qp] * 3.0 * deck.d_bar_E_phi_T * d_phi[j][qp] * d_phi[i][qp];

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
