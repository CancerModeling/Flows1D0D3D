////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

// Assembly class
void netfvfe::GradTafAssembly::assemble() {

  // Get required system alias
  // auto &grad_taf = d_model_p->get_grad_taf_assembly();
  auto &taf = d_model_p->get_taf_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();

  // Store current solution
  Gradient taf_grad;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    taf.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      taf_grad = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++)
        taf_grad.add_scaled(d_dphi[l][qp], taf.get_current_sol(l));

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        for (unsigned int d = 0; d < deck.d_dim; d++)
          d_Fe_var[d](i) += d_JxW[qp] * taf_grad(d) * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          for (unsigned int d = 0; d < deck.d_dim; d++)
            d_Ke_var[d][d](i, j) +=
              d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}
