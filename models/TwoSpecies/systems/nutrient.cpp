////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number twosp::initial_condition_nut(const Point &p, const Parameters &es,
                                      const std::string &system_name,
                                      const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void twosp::boundary_condition_nut(EquationSystems &es) {

  const auto *deck = es.parameters.get<InpDeck *>("input_deck");

  std::set<boundary_id_type> ids;
  if (deck->d_nutrient_bc_north)
    ids.insert(2);
  if (deck->d_nutrient_bc_south)
    ids.insert(0);
  if (deck->d_nutrient_bc_east)
    ids.insert(1);
  if (deck->d_nutrient_bc_west)
    ids.insert(3);

  // get variable ids
  auto &sys = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("nutrient"));

  // create constant function for bc and apply
  ConstFunction<Number> const_bc(1);
  DirichletBoundary diri_bc(ids, vars, &const_bc);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

// Assembly class
void twosp::NutAssembly::assemble() {
  assemble_1();
}

void twosp::NutAssembly::assemble_1() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;

  // Store current and old solution
  Real nut_old = 0.;
  Real tum_cur = 0.;

  Real tum_proj = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      nut_old = 0.;
      tum_cur = 0.;

      for (unsigned int l = 0; l < d_phi.size(); l++) {

        nut_old += d_phi[l][qp] * get_old_sol(l);
        tum_cur += d_phi[l][qp] * tum.get_current_sol_var(l, 0);
      }

      if (deck.d_assembly_method == 1) {

        compute_rhs =
            d_JxW[qp] * (nut_old + dt * deck.d_lambda_A * tum_cur);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_P * tum_cur);

      } else {

        tum_proj = util::project_concentration(tum_cur);

        compute_rhs =
            d_JxW[qp] * (nut_old + dt * deck.d_lambda_A * tum_proj);

        compute_mat = d_JxW[qp] * (1. + dt * deck.d_lambda_P * tum_proj);
      }

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe(i) += compute_rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke(i, j) += compute_mat * d_phi[j][qp] * d_phi[i][qp];

          // gradient term
          d_Ke(i, j) +=
              d_JxW[qp] * dt * deck.d_D_sigma * d_dphi[j][qp] * d_dphi[i][qp];
        }
      }

    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}