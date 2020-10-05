////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number avafv::initial_condition_hyp_kernel(const Point &p,
                                           const unsigned int &dim,
                                           const std::string &ic_type,
                                           const std::vector<double> &ic_center,
                                           const std::vector<double> &tum_ic_radius,
                                           const std::vector<double> &hyp_ic_radius) {

  const Point xc = util::to_point(ic_center);
  const Point dx = p - xc;

  if (ic_type == "tumor_hypoxic_spherical") {

    if (dx.norm() < tum_ic_radius[0] - 1.0E-12) {

      return 1. - util::exp_decay_function(
                    dx.norm() / tum_ic_radius[0], 4.);

    } else if (dx.norm() > tum_ic_radius[0] - 1.0E-12 and
               dx.norm() < hyp_ic_radius[0] - 1.0E-12) {

      return util::exp_decay_function(
        (dx.norm() - tum_ic_radius[0]) /
          (hyp_ic_radius[0] - tum_ic_radius[0]),
        4.);
    }
  } else if (ic_type == "tumor_hypoxic_elliptical") {

    // transform ellipse into ball of radius
    double small_ball_r = 0.;
    for (unsigned int i = 0; i < dim; i++)
      small_ball_r = tum_ic_radius[i] * tum_ic_radius[i];
    small_ball_r = std::sqrt(small_ball_r);

    Point p_small_ball =
      util::ellipse_to_ball(p, xc, tum_ic_radius, dim, small_ball_r);

    // transform ellipse into ball of radius
    double large_ball_r = 0.;
    for (unsigned int i = 0; i < dim; i++)
      large_ball_r = hyp_ic_radius[i] * hyp_ic_radius[i];
    large_ball_r = std::sqrt(large_ball_r);

    Point p_large_ball =
      util::ellipse_to_ball(p, xc, hyp_ic_radius, dim, large_ball_r);

    if (p_small_ball.norm() < small_ball_r - 1.0E-12) {

      return 1. -
             util::exp_decay_function(p_small_ball.norm() / small_ball_r, 4.);

    } else if (p_small_ball.norm() > small_ball_r - 1.0E-12 and
               p_small_ball.norm() < large_ball_r - 1.0E-12) {

      return util::exp_decay_function((p_small_ball.norm() - small_ball_r) /
                                        (large_ball_r - small_ball_r),
                                      4.);
    }
  }

  return 0.;
}

Number avafv::initial_condition_hyp(const Point &p, const Parameters &es,
                                    const std::string &system_name,
                                    const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Hypoxic");

  if (var_name == "hypoxic") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    double val = 0.;
    for (unsigned int i = 0; i < deck->d_tum_ic_data.size(); i++)
      val += initial_condition_hyp_kernel(p, deck->d_dim,
                                          deck->d_tum_ic_data[i].d_ic_type,
                                          deck->d_tum_ic_data[i].d_ic_center,
                                          deck->d_tum_ic_data[i].d_tum_ic_radius,
                                          deck->d_tum_ic_data[i].d_hyp_ic_radius);

    return val;
  }

  return 0.;
}

// Assembly class
void avafv::HypAssembly::assemble() {
  assemble_face();
  assemble_1();
}

void avafv::HypAssembly::assemble_face() {

  // Get required system alias
  // auto &hyp = d_model_p->get_hyp_assembly();
  auto &tum = d_model_p->get_tum_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Real> bc_vals;

  // to store pair of column dof and row-column matrix value
  std::vector<unsigned int> Ke_dof_row(1, 0);
  std::vector<Real> Ke_val_col;
  std::vector<unsigned int> Ke_dof_col;

  // local matrix and vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe(1);

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_hyp_neigh;
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real tum_cur = 0.;
  Real chem_tum_cur = 0.;
  Real hyp_proj = 0.;

  Real mobility_elem = 0.;

  // Store current and old solution of neighboring element
  Real tum_neigh_cur = 0.;
  Real chem_tum_neigh_cur = 0.;
  Real hyp_neigh_proj = 0.;

  Real mobility_neighbor = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution in this element
    tum_cur = tum.get_current_sol_var(0, 0);
    chem_tum_cur = tum.get_current_sol_var(0, 1);
    hyp_proj = util::project_concentration(get_current_sol(0));

    // face terms
    mobility_elem =
      deck.d_bar_M_H * pow(hyp_proj, 2) * pow(1. - hyp_proj, 2);

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);

        // get dof id
        // tum
        tum.init_var_dof(neighbor, dof_indices_tum_neigh,
                         dof_indices_tum_var_neigh);

        tum_neigh_cur =
          tum.get_current_sol_var(0, 0, dof_indices_tum_var_neigh);
        chem_tum_neigh_cur =
          tum.get_current_sol_var(0, 1, dof_indices_tum_var_neigh);

        // hyp
        init_dof(neighbor, dof_indices_hyp_neigh);
        hyp_neigh_proj = util::project_concentration(
          get_current_sol(0, dof_indices_hyp_neigh));

        // mobility in neighboring element
        mobility_neighbor = deck.d_bar_M_H *
                            pow(hyp_neigh_proj, 2) *
                            pow(1. - hyp_neigh_proj, 2);

        // mobility term due to div(grad(mu))
        // Goes to force
        compute_rhs = 0.;
        if (mobility_elem + mobility_neighbor > 1.0E-12)
          compute_rhs = dt * 2. * deck.d_face_by_h * mobility_elem *
                        mobility_neighbor /
                        (mobility_elem + mobility_neighbor);
        Fe(0) += -compute_rhs * (chem_tum_cur - chem_tum_neigh_cur);
      } // elem neighbor is not null
    }   // loop over faces

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}

void avafv::HypAssembly::assemble_1() {

  // Get required system alias
  // auto &hyp = d_model_p->get_hyp_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real hyp_old = 0.;
  Real hyp_cur = 0.;
  Real tum_cur = 0.;
  Real nut_cur = 0.;
  Real nec_cur = 0.;

  Real tum_proj = 0.;
  Real hyp_proj = 0.;
  Real nut_proj = 0.;
  Real nec_proj = 0.;

  Real mobility = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_dof(elem);
    nec.init_dof(elem);

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // get solution in this element
    nut_cur = nut.get_current_sol(0);
    hyp_old = get_old_sol(0);
    tum_cur = tum.get_current_sol_var(0, 0);
    nec_cur = nec.get_current_sol(0);

    // get projected solution
    hyp_proj = util::project_concentration(hyp_cur);

    mobility =
      deck.d_bar_M_H * pow(hyp_proj, 2) * pow(1. - hyp_proj, 2);

    if (deck.d_assembly_method == 1) {

      compute_rhs = deck.d_elem_size *
                    (hyp_old + dt * deck.d_lambda_PH * (tum_cur - nec_cur) *
                                 util::heaviside(deck.d_sigma_PH - nut_cur));

      compute_mat =
        deck.d_elem_size *
        (1. +
         dt *
           (deck.d_lambda_A +
            deck.d_lambda_PH * util::heaviside(deck.d_sigma_PH - nut_cur) +
            deck.d_lambda_HP * util::heaviside(nut_cur - deck.d_sigma_HP) +
            deck.d_lambda_HN * util::heaviside(deck.d_sigma_HN - nut_cur)));
    } else {

      tum_proj = util::project_concentration(tum_cur);
      nec_proj = util::project_concentration(nec_cur);
      nut_proj = util::project_concentration(nut_cur);

      compute_rhs = deck.d_elem_size *
                    (hyp_old + dt * deck.d_lambda_PH * (tum_proj - nec_proj) *
                                 util::heaviside(deck.d_sigma_PH - nut_proj));

      compute_mat =
        deck.d_elem_size *
        (1. +
         dt *
           (deck.d_lambda_A +
            deck.d_lambda_PH * util::heaviside(deck.d_sigma_PH - nut_proj) +
            deck.d_lambda_HP * util::heaviside(nut_proj - deck.d_sigma_HP) +
            deck.d_lambda_HN * util::heaviside(deck.d_sigma_HN - nut_proj)));
    }

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