////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"

void netfv::assemble_diffusion(BaseAssembly &sys) {

  // get tumor equation system
  EquationSystems &es = sys.d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Model parameters
  const auto &deck = sys.d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // get diffusion coefficient
  Real diff_coeff = 0.;
  if (sys.d_sys_name == "TAF")
    diff_coeff = deck.d_D_TAF;
  else if (sys.d_sys_name == "MDE")
    diff_coeff = deck.d_D_MDE;

  // to store pair of column dof and row-column matrix value
  std::vector<unsigned int> Ke_dof_row(1, 0);
  std::vector<Real> Ke_val_col;
  std::vector<unsigned int> Ke_dof_col;

  // local matrix and vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe(1);

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_sys_neigh;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    sys.init_dof(elem);

    // const unsigned int n_dofs = sys.d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = sys.get_global_dof_id(0);
    Fe(0) = 0.;

    // face terms
    {

      // loop over sides of the element
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);

          // get dof id
          // sys
          sys.init_dof(neighbor, dof_indices_sys_neigh);

          // diffusion
          const Real a_diff = dt * diff_coeff * deck.d_face_by_h;
          util::add_unique(sys.get_global_dof_id(0), a_diff, Ke_dof_col,
                           Ke_val_col);
          util::add_unique(dof_indices_sys_neigh[0], -a_diff, Ke_dof_col,
                           Ke_val_col);
        } // elem neighbor is not null
      }   // loop over faces

    } // terms over face of element

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    sys.d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    sys.d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}

void netfv::assemble_diffusion_advection(BaseAssembly &sys, PressureAssembly &pres, TumAssembly &tum) {

  // get tumor equation system
  EquationSystems &es = sys.d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Model parameters
  const auto &deck = sys.d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // get diffusion coefficient
  Real diff_coeff = 0.;
  if (sys.d_sys_name == "TAF")
    diff_coeff = deck.d_D_TAF;
  else if (sys.d_sys_name == "MDE")
    diff_coeff = deck.d_D_MDE;

  // to store pair of column dof and row-column matrix value
  std::vector<unsigned int> Ke_dof_row(1, 0);
  std::vector<Real> Ke_val_col;
  std::vector<unsigned int> Ke_dof_col;

  // local matrix and vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe(1);

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_sys_neigh;
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

    tum.init_var_dof(elem);
    pres.init_dof(elem);
    sys.init_dof(elem);

    // const unsigned int n_dofs = sys.d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = sys.get_global_dof_id(0);
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
          // sys
          sys.init_dof(neighbor, dof_indices_sys_neigh);

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
          const Real a_diff = dt * diff_coeff * deck.d_face_by_h;
          util::add_unique(sys.get_global_dof_id(0), a_diff, Ke_dof_col,
                           Ke_val_col);
          util::add_unique(dof_indices_sys_neigh[0], -a_diff, Ke_dof_col,
                           Ke_val_col);

          // advection
          Real mu_two_point = 0.;
          if (std::abs(chem_tum_cur + chem_tum_neigh_cur) > 1.E-12)
            mu_two_point = 2. * chem_tum_cur * chem_tum_neigh_cur /
                           (chem_tum_cur + chem_tum_neigh_cur);
          Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                     ((pres_cur - pres_neigh_cur) -
                      mu_two_point * (tum_cur - tum_neigh_cur));

          // upwinding
          if (v >= 0.)
            util::add_unique(sys.get_global_dof_id(0), dt * v, Ke_dof_col,
                             Ke_val_col);
          else
            util::add_unique(dof_indices_sys_neigh[0], dt * v, Ke_dof_col,
                             Ke_val_col);
        } // elem neighbor is not null
      }   // loop over faces

    } // terms over face of element

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    sys.d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    sys.d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}

void netfv::assemble_advection(BaseAssembly &sys, PressureAssembly &pres, TumAssembly &tum) {

  // get tumor equation system
  EquationSystems &es = sys.d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Model parameters
  const auto &deck = sys.d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // to store pair of column dof and row-column matrix value
  std::vector<unsigned int> Ke_dof_row(1, 0);
  std::vector<Real> Ke_val_col;
  std::vector<unsigned int> Ke_dof_col;

  // local matrix and vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe(1);

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_sys_neigh;
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

    tum.init_var_dof(elem);
    pres.init_dof(elem);
    sys.init_dof(elem);

    // const unsigned int n_dofs = sys.d_dof_indices_sys.size();

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = sys.get_global_dof_id(0);
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
          // sys
          sys.init_dof(neighbor, dof_indices_sys_neigh);

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

          // advection
          Real mu_two_point = 0.;
          if (std::abs(chem_tum_cur + chem_tum_neigh_cur) > 1.E-12)
            mu_two_point = 2. * chem_tum_cur * chem_tum_neigh_cur /
                           (chem_tum_cur + chem_tum_neigh_cur);
          Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                     ((pres_cur - pres_neigh_cur) -
                      mu_two_point * (tum_cur - tum_neigh_cur));

          // upwinding
          if (v >= 0.)
            util::add_unique(sys.get_global_dof_id(0), dt * v, Ke_dof_col,
                             Ke_val_col);
          else
            util::add_unique(dof_indices_sys_neigh[0], dt * v, Ke_dof_col,
                             Ke_val_col);
        } // elem neighbor is not null
      }   // loop over faces

    } // terms over face of element

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    sys.d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    sys.d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}
