////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_taf_source(const std::string &test_name, const Point &x, const
                      std::vector<double> &x0, const double &r) {

  if (test_name != "test_taf" and test_name != "test_taf_2")
    return 0.;

  double L_source_x = x0[0];
  double L_source_y = x0[1];
  double L_source_z = x0[2];
  const Point xc = Point(L_source_x, L_source_y, L_source_z);
  // spherical source
  if (false) {
    if ((x - xc).norm() < r)
      return 1.;
  } else {
    if (x(0) > L_source_x - r and x(0) < L_source_x + r and
        x(1) > L_source_y - r and x(1) < L_source_y + r) {

      //            compute_rhs += d_JxW[qp] * dt * deck->d_lambda_TAF *
      //                           std::sin(2. * (2. * M_PI) * x(2) / L_source);
      const Point x_plane = Point(x(0), x(1), 0.);
      const Point xc_plane = Point(L_source_x, L_source_y, 0.);
      if ((x_plane - xc_plane).norm() < r)
        return 1.;
    }
  }

  return 0.;
}
}

Number netfvfe::initial_condition_taf(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"TAF");

  return 0.;
}

// Assembly class
void netfvfe::TafAssembly::assemble() {
  assemble_1();
}

void netfvfe::TafAssembly::assemble_1() {

  // Get required system alias
  // auto &taf = d_model_p->get_ecm_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();  

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

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

    // init fe and element matrix and vector
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      taf_old = 0.; hyp_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        taf_old += d_phi[l][qp] * get_old_sol(l);
        hyp_cur += d_phi[l][qp] * hyp.get_current_sol(l);
      }

      if (deck.d_assembly_method == 1) {

        compute_rhs = d_JxW[qp] * (taf_old + dt * deck.d_lambda_TAF * hyp_cur);
        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_TAF * hyp_cur);

      } else {

        hyp_proj = util::project_concentration(hyp_cur);

        compute_rhs = d_JxW[qp] * (taf_old + dt * deck.d_lambda_TAF * hyp_proj);
        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_TAF * hyp_proj);
      }

      // add artificial source if any
      compute_rhs +=
          d_JxW[qp] * dt * deck.d_lambda_TAF *
          get_taf_source(deck.d_test_name, d_qpoints[qp],
                         deck.d_taf_source_center, deck.d_taf_source_radius);

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe(i) += compute_rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke(i, j) += compute_mat * d_phi[j][qp] * d_phi[i][qp];

          // gradient term
          d_Ke(i, j) +=
              d_JxW[qp] * dt * deck.d_D_TAF * d_dphi[j][qp] * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}