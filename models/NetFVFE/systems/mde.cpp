////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfvfe::initial_condition_mde(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"MDE");

  if (var_name == "mde") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    double val = 0.;
    for (unsigned int i=0; i<deck->d_tum_ic_data.size(); i++)
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
void netfvfe::MdeAssembly::assemble() {
  assemble_1();
}

void netfvfe::MdeAssembly::assemble_1() {

  // Get required system alias
  // auto &mde = d_model_p->get_mde_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

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

    // init fe and element matrix and vector
    init_fe(elem);

    // get finite-volume quantities
    nut_cur = nut.get_current_sol(0);
    nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      mde_old = 0.; tum_cur = 0.; nec_cur = 0.; ecm_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        mde_old += d_phi[l][qp] * get_old_sol(l);
        tum_cur += d_phi[l][qp] * tum.get_current_sol_var(l, 0);
        nec_cur += d_phi[l][qp] * nec.get_current_sol(l);
        ecm_cur += d_phi[l][qp] * ecm.get_current_sol(l);
      }

      if (deck.d_assembly_method == 1) {

        Real aux_1 = dt * deck.d_lambda_MDE_P * (tum_cur - nec_cur) * ecm_cur *
                     deck.d_sigma_HP / (1. + nut_cur);

        compute_rhs = d_JxW[qp] * (mde_old + aux_1);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_MDE_D + aux_1 +
                                   dt * deck.d_lambda_ECM_D * ecm_cur);
      } else {

        tum_proj = util::project_concentration(tum_cur);
        nec_proj = util::project_concentration(nec_cur);
        ecm_proj = util::project_concentration(ecm_cur);

        Real aux_1 = dt * deck.d_lambda_MDE_P * (tum_proj - nec_proj) *
                     ecm_proj * deck.d_sigma_HP / (1. + nut_proj);

        compute_rhs = d_JxW[qp] * (mde_old + aux_1);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_MDE_D + aux_1 +
                                   dt * deck.d_lambda_ECM_D * ecm_proj);
      }

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe(i) += compute_rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke(i, j) += compute_mat * d_phi[j][qp] * d_phi[i][qp];

          // gradient term
          d_Ke(i, j) +=
              d_JxW[qp] * dt * deck.d_D_MDE * d_dphi[j][qp] * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}