////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfv::initial_condition_tum(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"Tumor");

  if (var_name == "tumor") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    const unsigned int dim = deck->d_dim;
    const unsigned int num_ic = deck->d_tum_ic_data.size();
    if (num_ic == 0)
      return 0.;

    for (unsigned int ic = 0; ic < num_ic; ic++) {

      auto data = deck->d_tum_ic_data[ic];

      const std::string type =data.d_ic_type;
      const Point xc = Point(data.d_ic_center[0], data.d_ic_center[1],
                             data.d_ic_center[2]);
      const Point dx = p - xc;

      if (type == "tumor_spherical" or type == "tumor_hypoxic_spherical" or
          type == "tumor_spherical_sharp") {
        if (dx.norm() < data.d_tum_ic_radius[0] - 1.0E-12) {

          // out << "here tum ic\n";

          if (type == "tumor_spherical_sharp")
            return 1.;
          else
            return util::exp_decay_function(dx.norm() / data.d_tum_ic_radius[0],
                                          4.);
        }
      } else if (type == "tumor_elliptical" or
                 type == "tumor_hypoxic_elliptical") {

        // transform ellipse into ball of radius
        double small_ball_r = 0.;
        for (unsigned int i = 0; i < dim; i++)
          small_ball_r = data.d_tum_ic_radius[i] * data.d_tum_ic_radius[i];
        small_ball_r = std::sqrt(small_ball_r);

        Point p_small_ball = util::ellipse_to_ball(p, xc, data.d_tum_ic_radius,
                                                   dim, small_ball_r);

        if (p_small_ball.norm() < small_ball_r - 1.0E-12) {

          return util::exp_decay_function(p_small_ball.norm() / small_ball_r,
                                          4.);
        }
      }
    }

    return 0.;
  }

  return 0.;
}

// Assembly class
void netfv::TumAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else
    assemble_2();
}

void netfv::TumAssembly::assemble_face() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // to store pair of column dof and row-column matrix value
  std::vector<Real> Ke_phi_val_col;
  std::vector<unsigned int> Ke_phi_dof_col;
  std::vector<Real> Ke_mu_val_col;
  std::vector<unsigned int> Ke_mu_dof_col;

  std::vector<unsigned int> Ke_phi_dof_row(1);
  std::vector<unsigned int> Ke_mu_dof_row(1);

  // local matrix and vector
  DenseMatrix<Number> Ke_phi;
  DenseMatrix<Number> Ke_mu;

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_pres_neigh;
  std::vector<unsigned int> dof_indices_hyp_neigh;
  std::vector<unsigned int> dof_indices_nec_neigh;
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real pres_cur = 0.;
  Real tum_cur = 0.;
  Real chem_tum_cur = 0.;
  Real hyp_proj = 0.;
  Real pro_proj = 0.;

  // Store current and old solution of neighboring element
  Real pres_neigh_cur = 0.;
  Real tum_neigh_cur = 0.;
  Real chem_tum_neigh_cur = 0.;
  Real hyp_neigh_proj = 0.;
  Real pro_neigh_proj = 0.;

  Real mobility_elem = 0.;
  Real mobility_neighbor = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_var_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    pres.init_dof(elem);

    // const unsigned int n_dofs = tum.d_dof_indices_sys.size();
    // const unsigned int n_var_dofs = tum.d_dof_indices_sys_var[0].size();

    // reset matrix and force
    Ke_phi_dof_col.clear();
    Ke_phi_val_col.clear();
    Ke_mu_dof_col.clear();
    Ke_mu_val_col.clear();

    Ke_phi_dof_row[0] = get_var_global_dof_id(0, 0);
    Ke_mu_dof_row[0] = get_var_global_dof_id(0, 1);

    // get solution in this element
    pres_cur = pres.get_current_sol(0);
    tum_cur = get_current_sol_var(0, 0);
    chem_tum_cur = get_current_sol_var(0, 1);

    hyp_proj = util::project_concentration(hyp.get_current_sol(0));
    pro_proj = util::project_concentration(tum_cur - hyp.get_current_sol(0) - nec.get_current_sol(0));

    // face terms
    {
      mobility_elem =
          deck.d_bar_M_P * pow(pro_proj, 2) * pow(1. - pro_proj, 2) +
          deck.d_bar_M_H * pow(hyp_proj, 2) * pow(1. - hyp_proj, 2);

      // loop over sides of the element
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);

          // get dof id
          // pres
          pres.init_dof(neighbor, dof_indices_pres_neigh);
          pres_neigh_cur =
              pres.get_current_sol(0, dof_indices_pres_neigh);

          // tum
          init_var_dof(neighbor, dof_indices_tum_neigh, dof_indices_tum_var_neigh);

          tum_neigh_cur = get_current_sol_var(0, 0,
                                                       dof_indices_tum_var_neigh);
          chem_tum_neigh_cur = get_current_sol_var(0, 1,
                                                            dof_indices_tum_var_neigh);

          // hyp
          hyp.init_dof(neighbor, dof_indices_hyp_neigh);

          // nec
          nec.init_dof(neighbor, dof_indices_nec_neigh);

          // mobility in neighboring element
          hyp_neigh_proj = util::project_concentration(
              hyp.get_current_sol(0, dof_indices_hyp_neigh));
          pro_neigh_proj = util::project_concentration(
              tum_neigh_cur - hyp.get_current_sol(0, dof_indices_hyp_neigh) -
              nec.get_current_sol(0, dof_indices_nec_neigh));
          mobility_neighbor = deck.d_bar_M_P * pow(pro_neigh_proj, 2) *
                                  pow(1. - pro_neigh_proj, 2) +
                              deck.d_bar_M_H * pow(hyp_neigh_proj, 2) *
                                  pow(1. - hyp_neigh_proj, 2);

          // mobility term due to div(grad(mu))
          // Goes to row corresponding to phi and columns corresponding to mu
          Real a_mob = 0.;
          if (mobility_elem + mobility_neighbor > 1.E-12)
            a_mob = dt * 2. * deck.d_face_by_h * mobility_elem *
                             mobility_neighbor /
                             (mobility_elem + mobility_neighbor);
          util::add_unique(get_var_global_dof_id(0, 1), a_mob,
                           Ke_phi_dof_col, Ke_phi_val_col);
          util::add_unique(dof_indices_tum_var_neigh[1][0], -a_mob,
                           Ke_phi_dof_col, Ke_phi_val_col);

          // advection
          if (deck.d_advection_active) {
            // Goes to row and columns corresponding to phi
            Real mu_two_point = 0.;
            if (std::abs(chem_tum_cur + chem_tum_neigh_cur) > 1.0E-12)
              mu_two_point = 2. * chem_tum_cur * chem_tum_neigh_cur /
                             (chem_tum_cur + chem_tum_neigh_cur);
            Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                     ((pres_cur - pres_neigh_cur) -
                      mu_two_point * (tum_cur - tum_neigh_cur));

            // upwinding
            if (v >= 0.)
              util::add_unique(get_var_global_dof_id(0, 0), dt * v,
                               Ke_phi_dof_col, Ke_phi_val_col);
            else
              util::add_unique(dof_indices_tum_var_neigh[0][0], dt * v,
                               Ke_phi_dof_col, Ke_phi_val_col);
          }

          // interface term in cahn-hilliard
          // Goes to row corresponding to mu and columns corresponding to phi
          const Real a_inter = deck.d_epsilon_T * deck.d_epsilon_T *
              deck.d_face_by_h;
          util::add_unique(get_var_global_dof_id(0, 0), -a_inter,
                           Ke_mu_dof_col, Ke_mu_val_col);
          util::add_unique(dof_indices_tum_var_neigh[0][0], a_inter,
                           Ke_mu_dof_col, Ke_mu_val_col);
        } // elem neighbor is not null
      }   // loop over faces

    } // terms over face of element

    // add to matrix
    Ke_phi.resize(1, Ke_phi_dof_col.size());
    for (unsigned int i = 0; i < Ke_phi_dof_col.size(); i++)
      Ke_phi(0, i) = Ke_phi_val_col[i];

    d_sys.matrix->add_matrix(Ke_phi, Ke_phi_dof_row, Ke_phi_dof_col);

    Ke_mu.resize(1, Ke_mu_dof_col.size());
    for (unsigned int i = 0; i < Ke_mu_dof_col.size(); i++)
      Ke_mu(0, i) = Ke_mu_val_col[i];

    d_sys.matrix->add_matrix(Ke_mu, Ke_mu_dof_row, Ke_mu_dof_col);
  }
}

void netfv::TumAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Get required system alias
  // auto &tum = d_model_p->get_tum_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // to store dof of element
  std::vector<unsigned int> Ke_dof_row(2, 0);

  // local matrix and vector
  DenseMatrix<Number> Ke(2,2);
  DenseVector<Number> Fe(2);

  // Store current and old solution
  Real tum_old = 0.;
  Real tum_cur = 0.;
  Real nut_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;
  Real pro_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_var_dof(elem);
    nut.init_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);

    // const unsigned int n_dofs = tum.d_dof_indices_sys.size();
    // const unsigned int n_var_dofs = tum.d_dof_indices_sys_var[0].size();

    // reset matrix and force
    Ke_dof_row[0] = get_var_global_dof_id(0, 0);
    Ke_dof_row[1] = get_var_global_dof_id(0, 1);
    Ke.resize(2,2);
    Fe.resize(2);

    // get solution in this element
    nut_cur = nut.get_current_sol(0);
    tum_cur = get_current_sol_var(0, 0);
    hyp_cur = hyp.get_current_sol(0);
    nec_cur = nec.get_current_sol(0);
    pro_cur = tum_cur - hyp_cur - nec_cur;
    tum_old = get_old_sol_var(0, 0);

    // volume terms
    {
      //
      // phi term
      //
      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += tum_old * deck.d_elem_size;

      // source contribution to matrix
      Number a_source = deck.d_elem_size * dt * (deck.d_lambda_A + deck.d_lambda_P *
          nut_cur * pro_cur);
      Ke(0,0) += a_source;

      // source contribution to force
      Fe(0) +=
          deck.d_elem_size * dt *
          (deck.d_lambda_A * nec_cur + deck.d_lambda_P * nut_cur * pro_cur);

      //
      // mu term
      //
      // mass matrix
      Ke(1,1) += deck.d_elem_size;

      // implicit part of double-well
      // Goes to row corresponding to mu and column corresponding to phi
      Number a_dblw = -deck.d_elem_size * 3. * deck.d_bar_E_phi_T;
      Ke(1, 0) += a_dblw;

      // explicit part of double-well
      Fe(1) += deck.d_elem_size * deck.d_bar_E_phi_T * tum_old * (4.0 * pow(tum_old,
          2) - 6.0 * tum_old - 1.);

      // chemotactic term
      Fe(1) += -deck.d_elem_size * deck.d_chi_c * nut_cur;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void netfv::TumAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Get required system alias
  // auto &tum = d_model_p->get_tum_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // to store dof of element
  std::vector<unsigned int> Ke_dof_row(2, 0);

  // local matrix and vector
  DenseMatrix<Number> Ke(2,2);
  DenseVector<Number> Fe(2);

  // Store current and old solution
  Real tum_old = 0.;
  Real nut_proj = 0.;
  Real nec_proj = 0.;
  Real pro_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_var_dof(elem);
    nut.init_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);

    // const unsigned int n_dofs = tum.d_dof_indices_sys.size();
    // const unsigned int n_var_dofs = tum.d_dof_indices_sys_var[0].size();

    // reset matrix and force
    Ke_dof_row[0] = get_var_global_dof_id(0, 0);
    Ke_dof_row[1] = get_var_global_dof_id(0, 1);
    Ke.resize(2,2);
    Fe.resize(2);

    // get solution in this element
    tum_old = get_old_sol_var(0, 0);

    // get projected values of species
    pro_proj = util::project_concentration(get_current_sol_var(0, 0) -
        hyp.get_current_sol(0) - nec.get_current_sol(0));
    nut_proj = util::project_concentration(nut.get_current_sol(0));
    nec_proj = util::project_concentration(nec.get_current_sol(0));

    // volume terms
    {
      //
      // phi term
      //
      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += tum_old * deck.d_elem_size;

      // source contribution to matrix
      Number a_source = deck.d_elem_size * dt * (deck.d_lambda_A + deck.d_lambda_P *
                                                           nut_proj * pro_proj);
      Ke(0,0) += a_source;

      // source contribution to force
      Fe(0) +=
          deck.d_elem_size * dt *
          (deck.d_lambda_A * nec_proj + deck.d_lambda_P * nut_proj * pro_proj);

      //
      // mu term
      //
      // mass matrix
      Ke(1,1) += deck.d_elem_size;

      // implicit part of double-well
      // Goes to row corresponding to mu and column corresponding to phi
      Number a_dblw = -deck.d_elem_size * 3. * deck.d_bar_E_phi_T;
      Ke(1, 0) += a_dblw;

      // explicit part of double-well
      Fe(1) += deck.d_elem_size * deck.d_bar_E_phi_T * tum_old * (4.0 * pow(tum_old,
                                                                            2) - 6.0 * tum_old - 1.);

      // chemotactic term
      Fe(1) += -deck.d_elem_size * deck.d_chi_c * nut_proj;
    }

    // add to matrix
    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}