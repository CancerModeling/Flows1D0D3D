////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

void netfvfe::TumAssembly::solve_custom() {

  // Get required system alias
  auto &pro = d_model_p->get_pro_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  Real tum_cur = 0.;

  for (auto node : d_mesh.local_node_ptr_range()) {
    auto tum_dof = node->dof_number(d_sys.number(), 0, 0);
    auto pro_dof = node->dof_number(pro.d_sys.number(), 0, 0);
    auto hyp_dof = node->dof_number(hyp.d_sys.number(), 0, 0);
    auto nec_dof = node->dof_number(nec.d_sys.number(), 0, 0);

    tum_cur = pro.d_sys.current_solution(pro_dof) +
              hyp.d_sys.current_solution(hyp_dof) +
              nec.d_sys.current_solution(nec_dof);

    d_sys.solution->set(tum_dof, tum_cur);
  }

  d_sys.solution->close();
  *d_sys.old_local_solution = *d_sys.solution;
  d_sys.update();
}

// Assembly class
void netfvfe::TumAssembly::assemble() {
  assemble_1();
}

void netfvfe::TumAssembly::assemble_1() {

  // Get required system alias
  auto &pro = d_model_p->get_pro_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  // Store current and old solution
  Real pro_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;

  Real compute_rhs = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    pro.init_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      pro_cur = 0.;
      hyp_cur = 0.;
      nec_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        pro_cur += d_phi[l][qp] * pro.get_current_sol_var(l, 0);
        hyp_cur += d_phi[l][qp] * hyp.get_current_sol_var(l, 0);
        nec_cur += d_phi[l][qp] * nec.get_current_sol(l);
      }

      compute_rhs = d_JxW[qp] * (pro_cur + hyp_cur + nec_cur);

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe_var[0](i) += compute_rhs * d_phi[i][qp];

        d_Fe_var[1](i) += 0.;

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke_var[0][0](i, j) += d_phi[j][qp] * d_phi[i][qp];

          d_Ke_var[1][1](i, j) += d_phi[j][qp] * d_phi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(
      d_Ke, d_Fe, d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}