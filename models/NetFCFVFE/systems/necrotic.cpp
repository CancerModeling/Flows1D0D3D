////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfcfvfe::initial_condition_nec(const Point &p, const Parameters &es,
                                        const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Necrotic");

  return 0.;
}

// Assembly class
void netfcfvfe::NecAssembly::assemble() {
  assemble_1();
}

void netfcfvfe::NecAssembly::assemble_1() {

  // Get required system alias
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // Store current and old solution
  Real nec_old = 0.;
  Real nut_cur = 0.;
  Real hyp_cur = 0.;

  Real nut_proj = 0.;
  Real hyp_proj = 0.;

  Real compute_rhs = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    hyp.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    // get nutrient at this element
    nut_cur = nut.get_current_sol(0);
    nut_proj = util::project_concentration(nut_cur);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      nec_old = 0.;
      hyp_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        nec_old += d_phi[l][qp] * get_old_sol(l);
        hyp_cur += d_phi[l][qp] * hyp.get_current_sol_var(l, 0);
      }

      if (deck.d_assembly_method == 1) {
        compute_rhs =
          d_JxW[qp] *
          (nec_old + dt * deck.d_lambda_HN *
                       util::heaviside(deck.d_sigma_HN - nut_cur) *
                       hyp_cur);
      } else {

        hyp_proj = util::project_concentration(hyp_cur);

        compute_rhs =
          d_JxW[qp] *
          (nec_old + dt * deck.d_lambda_HN *
                       util::heaviside(deck.d_sigma_HN - nut_proj) *
                       hyp_proj);
      }

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe(i) += compute_rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke(i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}