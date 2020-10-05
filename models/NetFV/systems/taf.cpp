////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {

double get_taf_source(const std::string &test_name, const Point &x,
                      const std::vector<int> &type,
                      const std::vector<std::vector<double>> &centers,
                      const std::vector<double> &rads) {

  if (test_name != "test_taf" and test_name != "test_taf_2")
    return 0.;

  for (int i = 0; i < type.size(); i++) {

    const Point xc = util::to_point(centers[i]);
    auto d = (x - xc).norm();

    if (d < rads[i])
      return 1.;
  }

  return 0.;
}

} // namespace

Number netfv::initial_condition_taf(const Point &p, const Parameters &es,
                                    const std::string &system_name, const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "TAF");

  return 0.;
}

// Assembly class
void netfv::TafAssembly::assemble() {
  assemble_face();
  assemble_1();
}

void netfv::TafAssembly::assemble_face() {

  // call diffusion-advection calculation function
  if (d_model_p->get_input_deck().d_advection_active)
    netfv::assemble_diffusion_advection(*this,
                                        d_model_p->get_pres_assembly(),
                                        d_model_p->get_tum_assembly(), this->d_model_p);
  else
    netfv::assemble_diffusion(*this, this->d_model_p);
}

void netfv::TafAssembly::assemble_1() {

  // Get required system alias
  // auto &taf = d_model_p->get_ecm_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real taf_old = 0.;
  Real hyp_cur = 0.;

  Real hyp_proj = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    hyp.init_dof(elem);

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // get fields at this element
    hyp_cur = hyp.get_current_sol(0);
    taf_old = get_old_sol(0);

    if (deck.d_assembly_method == 1) {

      compute_rhs = deck.d_elem_size * (taf_old + dt * deck.d_lambda_TAF * hyp_cur);
      compute_mat = deck.d_elem_size * (1. + dt * deck.d_lambda_TAF * hyp_cur);

    } else {

      hyp_proj = util::project_concentration(hyp_cur);

      compute_rhs = deck.d_elem_size * (taf_old + dt * deck.d_lambda_TAF * hyp_proj);
      compute_mat = deck.d_elem_size * (1. + dt * deck.d_lambda_TAF * hyp_proj);
    }

    // add artificial source if any
    compute_rhs +=
      deck.d_elem_size * dt * deck.d_lambda_TAF *
      get_taf_source(deck.d_test_name, elem->centroid(),
                     deck.d_taf_source_type,
                     deck.d_taf_source_center, deck.d_taf_source_radius);

    // add
    Ke(0, 0) += compute_mat;
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