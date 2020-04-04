////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_nut_source(const std::string &test_name, const Point &x, const
std::vector<double> &x0, const double &r) {

  if (test_name != "test_tum_2")
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
    Point dx = Point(x(0), x(1), 0.) - xc;
    if (dx.norm() < r) {

      return 1.;
    }
  }

  return 0.;
}
}

Number netfv::initial_condition_nut(const Point &p, const Parameters &es,
                                     const std::string &system_name,
                                     const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<netfv::InputDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void netfv::boundary_condition_nut(EquationSystems &es) {

  const auto *deck = es.parameters.get<netfv::InputDeck *>("input_deck");

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
void netfv::NutAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_1d_coupling();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfv::NutAssembly::assemble_1d_coupling() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // coupling between tissue pressure and vessel pressure
  {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    DenseMatrix<Number> Ke(1,1);
    DenseVector<Number> Fe(1);

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;
      double c_v_k = pointer->c_v;

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = mesh.elem_ptr(e_id);
          pres.init_dof(elem);
          auto p_t_k = pres.get_current_sol(0);

          // get 3d nutrient
          init_dof(elem);
          auto c_t_k = get_current_sol(0);

          // implicit for c_t in source
          Ke(0, 0) = dt * factor_nut * deck.d_coupling_method_theta *
                      pointer->L_s[i] * J_b_data.half_cyl_surf * e_w;

          // explicit for c_v in source
          Fe(0) = dt * factor_nut * pointer->L_s[i] * J_b_data.half_cyl_surf *
                   e_w * (c_v_k - (1. - deck.d_coupling_method_theta) * c_t_k);

          // term due to pressure difference
          double c_transport = 0.;
          if (p_v_k - p_t_k >= 0.)
            c_transport = c_v_k;
          else
            c_transport = c_t_k;
          Fe(0) += dt * factor_nut * (1. - deck.d_osmotic_sigma) *
                   pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                   (p_v_k - p_t_k) * c_transport;

          // update matrix
          d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);
          d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
        }
      } // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }
}

void netfv::NutAssembly::assemble_face() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
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
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real pres_cur = 0.;
  Real tum_cur = 0.;
  Real chem_tum_cur = 0.;

  // Store current and old solution of neighboring element
  Real pres_neigh_cur = 0.;
  Real tum_neigh_cur = 0.;
  Real chem_tum_neigh_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_var_dof(elem);
    hyp.init_dof(elem);
    pres.init_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution in this element
    pres_cur = pres.get_current_sol(0);
    tum_cur = tum.get_current_sol_var(0, 0);
    chem_tum_cur = tum.get_current_sol_var(0, 1);

    // face terms
    {

      // loop over sides of the element
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);

          // get dof id
          // nut
          init_dof(neighbor, dof_indices_nut_neigh);

          // pres
          pres.init_dof(neighbor, dof_indices_pres_neigh);
          pres_neigh_cur =
              pres.get_current_sol(0, dof_indices_pres_neigh);

          // tum
          tum.init_var_dof(neighbor, dof_indices_tum_neigh, dof_indices_tum_var_neigh);
          tum_neigh_cur = tum.get_current_sol_var(0, 0,
                                                       dof_indices_tum_var_neigh);
          chem_tum_neigh_cur = tum.get_current_sol_var(0, 1,
                                                            dof_indices_tum_var_neigh);

          // diffusion
          const Real a_diff =
              factor_nut * dt * deck.d_D_sigma * deck.d_face_by_h;
          util::add_unique(get_global_dof_id(0), a_diff, Ke_dof_col,
                           Ke_val_col);
          util::add_unique(dof_indices_nut_neigh[0], -a_diff, Ke_dof_col,
                           Ke_val_col);

          // advection
          Real mu_two_point = 0.;
          if (std::abs(chem_tum_cur + chem_tum_neigh_cur) > 1.0E-12)
            mu_two_point = 2. * chem_tum_cur * chem_tum_neigh_cur /
                           (chem_tum_cur + chem_tum_neigh_cur);
          Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                     ((pres_cur - pres_neigh_cur) -
                      mu_two_point * (tum_cur - tum_neigh_cur));

          // upwinding
          if (v >= 0.)
            util::add_unique(get_global_dof_id(0), factor_nut * dt * v,
                             Ke_dof_col, Ke_val_col);
          else
            util::add_unique(dof_indices_nut_neigh[0], factor_nut * dt * v,
                             Ke_dof_col, Ke_val_col);

          // chemotactic term
          Number b_chi_c = factor_nut * dt * deck.d_chi_c * deck.d_D_sigma *
                           deck.d_face_by_h;
          Fe(0) += b_chi_c * (tum_cur - tum_neigh_cur);
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

void netfv::NutAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  //auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &taf = d_model_p->get_taf_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nut_old = 0.;
  Real tum_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;
  Real ecm_cur = 0.;
  Real mde_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_var_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    taf.init_dof(elem);
    ecm.init_dof(elem);
    mde.init_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get solution in this element
      nut_old = get_old_sol(0);
      tum_cur = tum.get_current_sol_var(0, 0);
      hyp_cur = hyp.get_current_sol(0);
      nec_cur = nec.get_current_sol(0);
      ecm_cur = ecm.get_current_sol(0);
      mde_cur = mde.get_current_sol(0);

      // mass matrix
      Ke(0,0) += factor_nut * deck.d_elem_size;

      // previous time step term
      Fe(0) += factor_nut * nut_old * deck.d_elem_size;

      // add source
      Fe(0) += factor_nut * deck.d_elem_size * dt *
               deck.d_lambda_ECM_D * ecm_cur * mde_cur;

      // handle all coupling with nutrient as source term
      Number a_source = deck.d_elem_size * dt *
                        (deck.d_lambda_P * (tum_cur - hyp_cur - nec_cur) +
                         deck.d_lambda_Ph * hyp_cur +
                         deck.d_lambda_ECM_P * (1. - ecm_cur) *
                             util::heaviside(ecm_cur - deck.d_bar_phi_ECM_P));
      Ke(0,0) += factor_nut * a_source;

      // add artificial source if asked
      double articial_source =
          get_nut_source(deck.d_test_name, elem->centroid(),
                         deck.d_nut_source_center, deck.d_nut_source_radius) -
          nut_old;
      if (articial_source > 0.)
        Fe(0) += deck.d_elem_size * dt * articial_source;
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

void netfv::NutAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

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
  const Real dt = es.parameters.get<Real>("time_step");
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nut_old = 0.;
  Real tum_proj = 0.;
  Real hyp_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;
  Real mde_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_var_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    taf.init_dof(elem);
    ecm.init_dof(elem);
    mde.init_dof(elem);

    const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;
    
    // volume terms
    {
      // get solution in this element
      nut_old = get_old_sol(0);

      // get projected values of species
      tum_proj = util::project_concentration(tum.get_current_sol_var(0, 0));
      hyp_proj = util::project_concentration(hyp.get_current_sol(0));
      nec_proj = util::project_concentration(nec.get_current_sol(0));
      ecm_proj = util::project_concentration(ecm.get_current_sol(0));
      mde_proj = util::project_concentration(mde.get_current_sol(0));

      // mass matrix
      Ke(0,0) += factor_nut * deck.d_elem_size;

      // previous time step term
      Fe(0) += factor_nut * nut_old * deck.d_elem_size;

      // add source
      Fe(0) += factor_nut * deck.d_elem_size * dt *
               deck.d_lambda_ECM_D * ecm_proj * mde_proj;

      // handle all coupling with nutrient as source term
      auto a_source = deck.d_elem_size * dt *
                        (deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) +
                         deck.d_lambda_Ph * hyp_proj +
                         deck.d_lambda_ECM_P * (1. - ecm_proj) *
                         util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));
      Ke(0,0) += factor_nut * a_source;

      // add artificial source if asked
      double articial_source = get_nut_source(deck.d_test_name, elem->centroid(),
                                              deck.d_nut_source_center, deck.d_nut_source_radius) - nut_old;
      if (articial_source > 0.)
        Fe(0) += deck.d_elem_size * dt * articial_source;
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

void netfv::NutAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  //auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &taf = d_model_p->get_taf_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nut_cur = 0.;
  Real nut_old = 0.;
  Real tum_proj = 0.;
  Real hyp_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;
  Real mde_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_var_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    taf.init_dof(elem);
    ecm.init_dof(elem);
    mde.init_dof(elem);

    // const unsigned int n_dofs = d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get solution in this element
      nut_cur = get_current_sol(0);
      nut_old = get_old_sol(0);

      // get projected values of species
      tum_proj = util::project_concentration(tum.get_current_sol_var(0, 0));
      hyp_proj = util::project_concentration(hyp.get_current_sol(0));
      nec_proj = util::project_concentration(nec.get_current_sol(0));
      ecm_proj = util::project_concentration(ecm.get_current_sol(0));
      mde_proj = util::project_concentration(mde.get_current_sol(0));

      // mass matrix
      Ke(0,0) += factor_nut * deck.d_elem_size;

      // previous time step term
      Fe(0) += factor_nut * nut_old * deck.d_elem_size;

      // add source
      Fe(0) += factor_nut * deck.d_elem_size * dt *
               deck.d_lambda_ECM_D * ecm_proj * mde_proj;

      // handle all coupling with nutrient as source term
      auto a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) +
                       deck.d_lambda_Ph * hyp_proj +
                       deck.d_lambda_ECM_P * (1. - ecm_proj) *
                       util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));

      // in this, instead of implicit, we consider it as a explicit forcing term
      Fe(0) += -factor_nut * a_source * nut_cur;

      // add artificial source if asked
      double articial_source = get_nut_source(deck.d_test_name, elem->centroid(),
                                              deck.d_nut_source_center, deck.d_nut_source_radius) - nut_old;
      if (articial_source > 0.)
        Fe(0) += deck.d_elem_size * dt * articial_source;
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