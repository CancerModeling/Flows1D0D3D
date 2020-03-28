////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number avafv::initial_condition_hyp_kernel(const Point &p,
                                            const avafv::InputDeck *deck) {

  const unsigned int dim = deck->d_dim;
  const unsigned int num_ic = deck->d_tum_ic_data.size();
  if (num_ic == 0)
    return 0.;

  for (unsigned int ic = 0; ic < num_ic; ic++) {

    auto data = deck->d_tum_ic_data[ic];

    const std::string type = data.d_ic_type;
    const Point xc =
        Point(data.d_ic_center[0], data.d_ic_center[1], data.d_ic_center[2]);
    const Point dx = p - xc;

    if (type == "tumor_hypoxic_spherical") {

      if (dx.norm() < data.d_tum_ic_radius[0] - 1.0E-12) {

        return 1. - util::exp_decay_function(
                        dx.norm() / data.d_tum_ic_radius[0], 4.);

      } else if (dx.norm() > data.d_tum_ic_radius[0] - 1.0E-12 and
                 dx.norm() < data.d_hyp_ic_radius[0] - 1.0E-12) {

        return util::exp_decay_function(
            (dx.norm() - data.d_tum_ic_radius[0]) /
                (data.d_hyp_ic_radius[0] - data.d_tum_ic_radius[0]),
            4.);
      }
    } else if (type == "tumor_hypoxic_elliptical") {

      // transform ellipse into ball of radius
      double small_ball_r = 0.;
      for (unsigned int i = 0; i < dim; i++)
        small_ball_r = data.d_tum_ic_radius[i] * data.d_tum_ic_radius[i];
      small_ball_r = std::sqrt(small_ball_r);

      Point p_small_ball =
          util::ellipse_to_ball(p, xc, data.d_tum_ic_radius, dim, small_ball_r);

      // transform ellipse into ball of radius
      double large_ball_r = 0.;
      for (unsigned int i = 0; i < dim; i++)
        large_ball_r = data.d_hyp_ic_radius[i] * data.d_hyp_ic_radius[i];
      large_ball_r = std::sqrt(large_ball_r);

      Point p_large_ball =
          util::ellipse_to_ball(p, xc, data.d_hyp_ic_radius, dim, large_ball_r);

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
  }

  return 0.;
}

Number avafv::initial_condition_hyp(const Point &p, const Parameters &es,
                                     const std::string &system_name,
                                     const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Hypoxic");

  if (var_name == "hypoxic") {

    const auto *deck = es.get<avafv::InputDeck *>("input_deck");

    return initial_condition_hyp_kernel(p, deck);
  }

  return 0.;
}

// Assembly class
void avafv::HypAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  assemble_vol();
}

void avafv::HypAssembly::assemble_face() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &hyp = d_model_p->get_hyp_assembly();
  auto &tum = d_model_p->get_tum_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

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
  std::vector<unsigned int> dof_indices_pres_neigh;
  std::vector<unsigned int> dof_indices_hyp_neigh;
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real chem_tum_cur = 0.;
  Real hyp_proj = 0.;

  // Store current and old solution of neighboring element
  Real chem_tum_neigh_cur = 0.;
  Real hyp_neigh_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_var_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution in this element
    chem_tum_cur = tum.get_current_sol_var(0, 1);
    hyp_proj = util::project_concentration(get_current_sol(0));

    // face terms
    {
      const Real mobility_elem =
          deck.d_bar_M_H * pow(hyp_proj, 2) * pow(1. - hyp_proj, 2);

      // loop over sides of the element
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);

          // get dof id
          // tum
          tum.init_var_dof(neighbor, dof_indices_tum_neigh,
                           dof_indices_tum_var_neigh);

          chem_tum_neigh_cur =
              tum.get_current_sol_var(0, 1, dof_indices_tum_var_neigh);

          // hyp
          init_dof(neighbor, dof_indices_hyp_neigh);
          hyp_neigh_proj = util::project_concentration(
              get_current_sol(0, dof_indices_hyp_neigh));

          // mobility in neighboring element
          const Real mobility_neighbor = deck.d_bar_M_H *
                                         pow(hyp_neigh_proj, 2) *
                                         pow(1. - hyp_neigh_proj, 2);

          // mobility term due to div(grad(mu))
          // Goes to force
          Real a_mob = 0.;
          if (mobility_elem + mobility_neighbor > 1.0E-12)
            a_mob = dt * 2. * deck.d_face_by_h * mobility_elem *
                               mobility_neighbor /
                               (mobility_elem + mobility_neighbor);
          Fe(0) += -a_mob * (chem_tum_cur - chem_tum_neigh_cur);
        } // elem neighbor is not null
      }   // loop over faces

    } // terms over face of element

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}

void avafv::HypAssembly::assemble_vol() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &hyp = d_model_p->get_hyp_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real hyp_old = 0.;
  Real nut_proj = 0.;
  Real via_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    tum.init_var_dof(elem);
    nec.init_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get solution in this element
      hyp_old = get_old_sol(0);

      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      via_proj = util::project_concentration(tum.get_current_sol_var(0, 0) -
                                             nec.get_current_sol(0));

      // mass matrix
      Ke(0, 0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += hyp_old * deck.d_elem_size;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_PH * via_proj *
               util::heaviside(deck.d_sigma_PH - nut_proj);

      // matrix contribution from source
      Number a_source =
          deck.d_elem_size * dt *
          (deck.d_lambda_A +
           deck.d_lambda_PH * util::heaviside(deck.d_sigma_PH - nut_proj) +
           deck.d_lambda_HP * util::heaviside(nut_proj - deck.d_sigma_HP) +
           deck.d_lambda_HN * util::heaviside(deck.d_sigma_HN - nut_proj));
      Ke(0, 0) += a_source;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}