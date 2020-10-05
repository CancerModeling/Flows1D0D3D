////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number avafv::initial_condition_nec(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Necrotic");

  return 0.;
}

// Assembly class
void avafv::NecAssembly::assemble() {
  assemble_1();
}

void avafv::NecAssembly::assemble_1() {

  // Get required system alias
  // auto &nec = d_model_p->get_nec_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

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

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // get solution in this element
    nut_cur = nut.get_current_sol(0);
    nec_old = get_old_sol(0);
    hyp_cur = hyp.get_current_sol(0);

    if (deck.d_assembly_method == 1) {
      compute_rhs =
        deck.d_elem_size *
        (nec_old + dt * deck.d_lambda_HN *
                     util::heaviside(deck.d_sigma_HN - nut_cur) *
                     hyp_cur);
    } else {

      hyp_proj = util::project_concentration(hyp_cur);
      nut_proj = util::project_concentration(nut_cur);

      compute_rhs =
        deck.d_elem_size *
        (nec_old + dt * deck.d_lambda_HN *
                     util::heaviside(deck.d_sigma_HN - nut_proj) *
                     hyp_proj);
    }

    // add
    Ke(0, 0) += deck.d_elem_size;

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