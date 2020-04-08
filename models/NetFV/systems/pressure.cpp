////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {

double get_exact_source(const Point &p) {

  return std::sin(2. * M_PI * p(0)) * std::sin(2. * M_PI * p(1)) *
         std::sin(2.* M_PI * p(2));
}
}

Number netfv::initial_condition_pres(const Point &p, const Parameters &es,
                                      const std::string &system_name,
                                      const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Pressure");

  if (var_name == "pressure") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_pressure_ic_val;
  }

  return 0.;
}

void netfv::boundary_condition_pres(EquationSystems &es) {

  const auto *deck = es.parameters.get<InpDeck *>("input_deck");

  std::set<boundary_id_type> ids;
  if (deck->d_nutrient_bc_north)
    ids.insert(2);
  if (deck->d_nutrient_bc_north)
    ids.insert(0);
  if (deck->d_nutrient_bc_north)
    ids.insert(1);
  if (deck->d_nutrient_bc_north)
    ids.insert(3);

  // get variable ids
  auto &sys = es.get_system<TransientLinearImplicitSystem>("Pressure");
  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("pressure"));

  // create constant function for bc and apply
  ConstFunction<Number> cf(deck->d_pressure_bc_val);
  DirichletBoundary diri_bc(ids, vars, &cf);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

// Assembly class

void netfv::PressureAssembly::assemble() {
  assemble_1d_coupling();
  assemble_face();
  assemble_1();
}

void netfv::PressureAssembly::assemble_1d_coupling() {

  // Get required system alias
  // auto &pres = d_model_p->get_pres_assembly();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();
  const double factor_p = deck.d_assembly_factor_p_t;

  // coupling between tissue pressure and vessel pressure
  {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    DenseMatrix<Number> Ke(1,1);
    DenseVector<Number> Fe(1);

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = d_mesh.elem_ptr(e_id);
          init_dof(elem);
          auto p_t_k = get_current_sol(0);

          // implicit for p_t in source
          Ke(0,0) = factor_p * deck.d_coupling_method_theta * pointer->L_p[i] *
                     J_b_data.half_cyl_surf *
                     e_w;

          // explicit for p_v in source
          Fe(0) = factor_p * pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                   (p_v_k - (1. - deck.d_coupling_method_theta) * p_t_k);

          //          if (pointer->index < 10 and i == 0 and e < 5)
          //            out << "index: " << pointer->index << ", p_v: " << p_v_k
          //                << ", p_t: " << p_t_k << ", Ke: " << Ke(0, 0)
          //                << ", Fe: " << Fe(0) << "\n";

          // update matrix
          d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);
          d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
        }
      } // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }
}

void netfv::PressureAssembly::assemble_face() {

  // Get required system alias
  // auto &pres = d_model_p->get_pres_assembly();
  auto &tum = d_model_p->get_tum_assembly();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real factor_p = deck.d_assembly_factor_p_t;

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
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real tum_cur = 0.;
  Real chem_tum_cur = 0.;

  // Store current and old solution of neighboring element
  Real tum_neigh_cur = 0.;
  Real chem_tum_neigh_cur = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto & elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);

    // reset matrix anf force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution at this element
    tum_cur = tum.get_current_sol_var(0, 0);
    chem_tum_cur = tum.get_current_sol_var(0, 1);

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);

        // get dof id
        // pres
        init_dof(neighbor, dof_indices_pres_neigh);

        // tum
        tum.init_var_dof(neighbor, dof_indices_tum_neigh,
                         dof_indices_tum_var_neigh);

        // get coefficient
        compute_mat = factor_p * deck.d_tissue_flow_rho *
                            deck.d_tissue_flow_coeff * deck.d_face_by_h;

        // diffusion
        util::add_unique(get_global_dof_id(0), compute_mat, Ke_dof_col,
                         Ke_val_col);
        util::add_unique(dof_indices_pres_neigh[0], -compute_mat, Ke_dof_col,
                         Ke_val_col);

        // div(chem_tum * grad(tum)) term
        tum_neigh_cur =
            tum.get_current_sol_var(0, 0, dof_indices_tum_var_neigh);
        chem_tum_neigh_cur =
            tum.get_current_sol_var(0, 1, dof_indices_tum_var_neigh);

        compute_rhs = 0.;
        if (std::abs(chem_tum_cur + chem_tum_neigh_cur) > 1.0e-12)
          compute_rhs = 2. * compute_mat * chem_tum_cur * chem_tum_neigh_cur /
                             (chem_tum_cur + chem_tum_neigh_cur);

        Fe(0) += compute_rhs * (tum_cur - tum_neigh_cur);
      } // elem neighbor is not null
    }   // loop over faces

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  } // element loop
}

void netfv::PressureAssembly::assemble_1() {

  // Get required system alias
  // auto &pres = d_model_p->get_pres_assembly();

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}