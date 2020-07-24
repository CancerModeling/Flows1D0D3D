////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_nut_source(const std::string &test_name, const Point &x,
                      const std::vector<double> &x0, const double &r) {

  if (test_name != "test_tum_2")
    return 0.;

  // For cylinder, take z-coord as 0
  const Point xc = util::to_point({x0[0], x0[1], 0.});

  // spherical source
  Point dx = Point(x(0), x(1), 0.) - xc;
  if (dx.norm() < r)
    return 1.;

  return 0.;
}
} // namespace

Number netfvfe::initial_condition_nut(const Point &p, const Parameters &es,
                                      const std::string &system_name,
                                      const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void netfvfe::boundary_condition_nut(EquationSystems &es) {

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
  ConstFunction<Number> one(1);
  DirichletBoundary diri_bc(ids, vars, &one);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

// Assembly class
void netfvfe::NutAssembly::assemble() {
  assemble_1d_coupling();
  assemble_face();
  assemble_1();
}

void netfvfe::NutAssembly::assemble_1d_coupling() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // coupling between tissue pressure and vessel pressure
  {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    DenseMatrix<Number> Ke(1, 1);
    DenseVector<Number> Fe(1);

    double h_3D = network.h_3D;
    int N_3D = network.N_3D;

    while (pointer) {

      std::vector<double> coord = pointer->coord;

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;
      double c_v_k = pointer->c_v;

      for (int i = 0; i < numberOfNeighbors; i++) {

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord, coord_neighbor);
        double L_s = pointer->L_s[i];
        double L_p = pointer->L_p[i];

        // Coupling terms
        int N_s = deck.d_num_points_length;
        int N_theta = deck.d_num_points_angle;

        std::vector<double> weights;
        std::vector<int> id_3D_elements;

        double length_edge = 0.5 * length;

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        util::unet::determineWeightsAndIds(
            N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D,
            length_edge, weights, id_3D_elements, d_mesh, true);

        // Add coupling entry
        int numberOfElements = id_3D_elements.size();

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {

            // get 3d pressure
            const auto *elem = d_mesh.elem_ptr(id_3D_elements[j]);
            if (elem->processor_id() == d_model_p->get_comm()->rank()) {
              pres.init_dof(elem);
              auto p_t_k = pres.get_current_sol(0);

              // get 3d nutrient
              init_dof(elem);
              auto c_t_k = get_current_sol(0);

              // implicit for c_t in source
              Ke(0, 0) = dt * factor_nut * L_s * surface_area * weights[j];

              // explicit for c_v in source
              Fe(0) = dt * factor_nut * L_s * surface_area * weights[j] * c_v_k;

              // osmotic reflection term
              if (p_v_k - p_t_k > 0.0) {

                // 3D equation
                // 2pi R (p_v - p_t) phi_v term in right hand side of 3D equation
                Fe(0) += dt * factor_nut * (1. - deck.d_osmotic_sigma) * L_p *
                         surface_area * weights[j] * (p_v_k - p_t_k) * c_v_k;

              } else {

                // 3D equation
                // 2pi R (p_v - p_t) phi_sigma term in right hand side of 3D equation
                Ke(0, 0) += -dt * factor_nut * (1. - deck.d_osmotic_sigma) *
                            L_p * surface_area * weights[j] * (p_v_k - p_t_k);
              }

              // update matrix
              d_sys.matrix->add_matrix(Ke, d_dof_indices_sys,
                                       d_dof_indices_sys);
              d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
            }
          }
        } // loop over 3D elements
      }   // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }
}

void netfvfe::NutAssembly::assemble_face() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

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
  std::vector<unsigned int> dof_indices_nut_neigh;
  std::vector<unsigned int> dof_indices_pres_neigh;

  // Store current and old solution
  Real pres_cur = 0.;
  Real chem_tum_cur = 0.;

  Gradient tum_grad = 0.;

  // Store current and old solution of neighboring element
  Real pres_neigh_cur = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);
    pres.init_dof(elem);

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution in this element
    pres_cur = pres.get_current_sol(0);

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);

        // get dof id
        // nut
        init_dof(neighbor, dof_indices_nut_neigh);

        // pres
        pres.init_dof(neighbor, dof_indices_pres_neigh);
        pres_neigh_cur = pres.get_current_sol(0, dof_indices_pres_neigh);

        // diffusion
        const Real a_diff = factor_nut * dt * deck.d_D_sigma * deck.d_face_by_h;
        util::add_unique(get_global_dof_id(0), a_diff, Ke_dof_col, Ke_val_col);
        util::add_unique(dof_indices_nut_neigh[0], -a_diff, Ke_dof_col,
                         Ke_val_col);

        // advection (pressure gradient term)
        Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                 (pres_cur - pres_neigh_cur);

        // upwinding
        if (v >= 0.)
          util::add_unique(get_global_dof_id(0), factor_nut * dt * v,
                           Ke_dof_col, Ke_val_col);
        else
          util::add_unique(dof_indices_nut_neigh[0], factor_nut * dt * v,
                           Ke_dof_col, Ke_val_col);

        // advection (mu_T grad(phi_T) term) and chemotactic term
        // these terms require integration over face of an element
        tum.d_fe_face->reinit(elem, side);

        // loop over quadrature points
        for (unsigned int qp = 0; qp < tum.d_qrule_face.n_points(); qp++) {

          chem_tum_cur = 0.;
          tum_grad = 0.;
          for (unsigned int l = 0; l < tum.d_phi_face.size(); l++) {

            chem_tum_cur +=
                tum.d_phi_face[l][qp] * tum.get_current_sol_var(l, 1);

            tum_grad.add_scaled(tum.d_dphi_face[l][qp],
                                tum.get_current_sol_var(l, 0));
          }

          // chemotactic term
          Fe(0) += -factor_nut * tum.d_JxW_face[qp] * dt * deck.d_D_sigma *
                   deck.d_chi_c * tum_grad * tum.d_qface_normals[qp];

          // advection term
          Real v_mu = factor_nut * tum.d_JxW_face[qp] * dt *
                      deck.d_tissue_flow_coeff * chem_tum_cur *
                      (tum_grad * tum.d_qface_normals[qp]);

          // goes to the dof of element (not the neighbor)
          util::add_unique(get_global_dof_id(0), v_mu, Ke_dof_col, Ke_val_col);
        } // loop over quadrature points on face

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

void netfvfe::NutAssembly::assemble_1() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &taf = d_model_p->get_taf_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // Store current and old solution
  Real nut_old = 0.;
  Real tum_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;
  Real ecm_cur = 0.;
  Real mde_cur = 0.;

  Real tum_proj = 0.;
  Real hyp_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;
  Real mde_proj = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    taf.init_dof(elem);
    ecm.init_dof(elem);
    mde.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);
    hyp.init_fe(elem);

    // get finite-volume quantities
    nut_old = get_old_sol(0);

    // add finite-volume contribution to matrix and vector
    d_Fe(0) += deck.d_elem_size * factor_nut * nut_old;
    d_Ke(0, 0) += deck.d_elem_size * factor_nut;

    // for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
    for (unsigned int qp = 0; qp < hyp.d_qrule.n_points(); qp++) {

      // Computing solution
      tum_cur = 0.;
      hyp_cur = 0.;
      nec_cur = 0.;
      ecm_cur = 0.;
      mde_cur = 0.;

      // for (unsigned int l = 0; l < d_phi.size(); l++) {
      for (unsigned int l = 0; l < hyp.d_phi.size(); l++) {

        tum_cur += hyp.d_phi[l][qp] * tum.get_current_sol_var(l, 0);
        hyp_cur += hyp.d_phi[l][qp] * hyp.get_current_sol(l);
        nec_cur += hyp.d_phi[l][qp] * nec.get_current_sol(l);
        ecm_cur += hyp.d_phi[l][qp] * ecm.get_current_sol(l);
        mde_cur += hyp.d_phi[l][qp] * mde.get_current_sol(l);
      }

      if (deck.d_assembly_method == 1) {

        compute_rhs =
            hyp.d_JxW[qp] * dt * deck.d_lambda_ECM_D * ecm_cur * mde_cur;

        compute_mat = hyp.d_JxW[qp] * dt *
                      (deck.d_lambda_P * (tum_cur - hyp_cur - nec_cur) +
                       deck.d_lambda_Ph * hyp_cur +
                       deck.d_lambda_ECM_P * (1. - ecm_cur) *
                           util::heaviside(ecm_cur - deck.d_bar_phi_ECM_P));

      } else {

        mde_proj = util::project_concentration(mde_cur);
        ecm_proj = util::project_concentration(ecm_cur);
        tum_proj = util::project_concentration(tum_cur);
        hyp_proj = util::project_concentration(hyp_cur);
        nec_proj = util::project_concentration(nec_cur);

        compute_rhs =
            hyp.d_JxW[qp] * dt * deck.d_lambda_ECM_D * ecm_proj * mde_proj;

        compute_mat = hyp.d_JxW[qp] * dt *
                      (deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) +
                       deck.d_lambda_Ph * hyp_proj +
                       deck.d_lambda_ECM_P * (1. - ecm_proj) *
                           util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));
      }

      // add artificial source if asked
      Real artificial_source =
          get_nut_source(deck.d_test_name, d_qpoints[qp],
                         deck.d_nut_source_center, deck.d_nut_source_radius) -
          nut_old;
      if (artificial_source > 0.)
        compute_rhs += hyp.d_JxW[qp] * dt * artificial_source;

      // add rhs
      d_Fe(0) += factor_nut * compute_rhs;

      // add matrix
      d_Ke(0, 0) += factor_nut * compute_mat;

    } // loop over quadrature points

    // add to matrix
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}