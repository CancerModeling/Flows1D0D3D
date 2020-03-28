////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "systems.hpp"
#include "../model.hpp"

void netfc::VelocityAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  if (deck.d_dim > 2)
    assemble_3d();
  else
    assemble_2d();
}

void netfc::VelocityAssembly::assemble_3d() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // Velocity
  auto &vel =
      es.get_system<TransientLinearImplicitSystem>("Velocity");
  std::vector<unsigned int> v_vel(3);
  v_vel[0] = vel.variable_number("velocity_x");
  v_vel[1] = vel.variable_number("velocity_y");
  v_vel[2] = vel.variable_number("velocity_z");

  const DofMap &vel_map = vel.get_dof_map();
  std::vector<unsigned int> dof_indices_vel;
  std::vector<std::vector<dof_id_type>> dof_indices_vel_var(3);

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type_vel = vel.variable_type(0);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_type_vel));
  QGauss qrule(dim, fe_type_vel.default_quadrature_order());
  fe_vel->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real hyd_cond = deck.d_tissue_flow_K / deck.d_tissue_flow_mu;
  const double factor_p_t = deck.d_coupling_factor_p_t;

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[3][3] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[3] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  std::vector<dof_id_type> dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    vel_map.dof_indices(elem, dof_indices_vel);
    for (unsigned int var = 0; var < 3; var++)
      vel_map.dof_indices(elem, dof_indices_vel_var[var], v_vel[var]);

    pres_map.dof_indices(elem, dof_indices_pres, v_pres);

    const unsigned int n_dofs = dof_indices_vel.size();
    const unsigned int n_var_dofs = dof_indices_vel_var[0].size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 3; var_i++)
      for (unsigned int var_j = 0; var_j < 3; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 3; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    Real p_cur = pres.current_solution(dof_indices_pres[0]);
    Gradient grad_p_cur;
    // Compute grad(p), which is constant in the element, at the center of
    // the element
    {
      unsigned int nx = 0, ny = 0, nz = 0;
      // loop over sides of the element
      const auto elem_center = elem->centroid();
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);
          const auto neighbor_center = neighbor->centroid();
          const auto dx = elem_center - neighbor_center;

          // get dof id
          std::vector<unsigned int> dof_indices_pres_neigh;
          pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);
          Real p_neigh_cur = pres.current_solution(dof_indices_pres_neigh[0]);

          // check if dx is aligned along x or y or z direction
          std::string dx_align = util::get_vec_major_axis(dx);

          if (dx_align == "x") {
            nx++;
            grad_p_cur(0) += (p_cur - p_neigh_cur) / dx(0);
          } else if (dx_align == "y") {
            ny++;
            grad_p_cur(1) += (p_cur - p_neigh_cur) / dx(1);
          } else if (dx_align == "z") {
            nz++;
            grad_p_cur(2) += (p_cur - p_neigh_cur) / dx(2);
          }
        }
      } // loop over neighbors

      // divide by number of neighbors
      grad_p_cur(0) = grad_p_cur(0) / nx;
      grad_p_cur(1) = grad_p_cur(1) / ny;
      grad_p_cur(2) = grad_p_cur(2) / nz;
    }

    // loop over quadrature points
    {
      for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

        Number chem_tum_cur = 0.;
        Gradient tum_grad;

        for (unsigned int l = 0; l < phi.size(); l++) {

          chem_tum_cur +=
              phi[l][qp] * tum.current_solution(dof_indices_tum_var[1][l]);
          tum_grad.add_scaled(dphi[l][qp],
                              tum.current_solution(dof_indices_tum_var[0][l]));
        }

        // Debug
        chem_tum_cur = 0.;
        tum_grad = 0.;

        Gradient compute_rhs =
            JxW[qp] * factor_p_t *
            (-hyd_cond * (grad_p_cur - chem_tum_cur * tum_grad));

        // Assembling matrix
        for (unsigned int i = 0; i < phi.size(); i++) {

          Fe_var[0](i) += compute_rhs(0) * phi[i][qp];
          Fe_var[1](i) += compute_rhs(1) * phi[i][qp];
          Fe_var[2](i) += compute_rhs(2) * phi[i][qp];

          for (unsigned int j = 0; j < phi.size(); j++) {

            Ke_var[0][0](i, j) +=
                JxW[qp] * factor_p_t * phi[j][qp] * phi[i][qp];
            Ke_var[1][1](i, j) +=
                JxW[qp] * factor_p_t * phi[j][qp] * phi[i][qp];
            Ke_var[2][2](i, j) +=
                JxW[qp] * factor_p_t * phi[j][qp] * phi[i][qp];
          }
        }
      } // loop over quadrature points

      vel_map.heterogenously_constrain_element_matrix_and_vector(
          Ke, Fe, dof_indices_vel);
      vel.matrix->add_matrix(Ke, dof_indices_vel);
      vel.rhs->add_vector(Fe, dof_indices_vel);
    }
  }
}

void netfc::VelocityAssembly::assemble_2d() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // Velocity
  auto &vel =
      es.get_system<TransientLinearImplicitSystem>("Velocity");
  std::vector<unsigned int> v_vel(2);
  v_vel[0] = vel.variable_number("velocity_x");
  v_vel[1] = vel.variable_number("velocity_y");

  const DofMap &vel_map = vel.get_dof_map();
  std::vector<unsigned int> dof_indices_vel;
  std::vector<std::vector<dof_id_type>> dof_indices_vel_var(2);

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type_vel = vel.variable_type(0);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_type_vel));
  QGauss qrule(dim, fe_type_vel.default_quadrature_order());
  fe_vel->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real hyd_cond = deck.d_tissue_flow_K / deck.d_tissue_flow_mu;

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[2][2] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] = {DenseSubVector<Number>(Fe),
                                      DenseSubVector<Number>(Fe)};

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    vel_map.dof_indices(elem, dof_indices_vel);
    for (unsigned int var = 0; var < 3; var++)
      vel_map.dof_indices(elem, dof_indices_vel_var[var], v_vel[var]);

    pres_map.dof_indices(elem, dof_indices_pres, v_pres);

    const unsigned int n_dofs = dof_indices_vel.size();
    const unsigned int n_var_dofs = dof_indices_vel_var[0].size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2; var_i++)
      for (unsigned int var_j = 0; var_j < 2; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);

    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    Real p_cur = pres.current_solution(dof_indices_pres[0]);
    Gradient grad_p_cur;
    // Compute grad(p), which is constant in the element, at the center of
    // the element
    {
      unsigned int nx = 0, ny = 0;
      // loop over sides of the element
      const auto elem_center = elem->centroid();
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem * neighbor = elem->neighbor_ptr(side);
          const auto neighbor_center = neighbor->centroid();
          const auto dx = elem_center - neighbor_center;

          // get dof id
          std::vector<unsigned int> dof_indices_pres_neigh;
          pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);
          Real p_neigh_cur = pres.current_solution(dof_indices_pres_neigh[0]);

          // check if dx is aligned along x or y or z direction
          std::string dx_align = util::get_vec_major_axis(dx);

          if (dx_align == "x") {
            nx++;
            grad_p_cur(0) += (p_cur - p_neigh_cur) / dx.norm();
          } else if (dx_align == "y") {
            ny++;
            grad_p_cur(1) += (p_cur - p_neigh_cur) / dx.norm();
          }
        }
      } // loop over neighbors

      // normalize gradient
      grad_p_cur(0) = grad_p_cur(0) / nx;
      grad_p_cur(1) = grad_p_cur(1) / ny;
    }

    // loop over quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      Number chem_tum_cur = 0.;
      Gradient tum_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {

        chem_tum_cur += phi[l][qp] *
                        tum.current_solution(dof_indices_tum_var[1][l]);
        tum_grad.add_scaled(dphi[l][qp],
                            tum.current_solution(dof_indices_tum_var[0][l]));
      }

      Gradient compute_rhs =
          JxW[qp] * (-hyd_cond * (grad_p_cur - chem_tum_cur * tum_grad));

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fe_var[0](i) += compute_rhs(0) * phi[i][qp];
        Fe_var[1](i) += compute_rhs(1) * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke_var[0][0](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
          Ke_var[1][1](i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    vel_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe,
                                                               dof_indices_vel);
    vel.matrix->add_matrix(Ke, dof_indices_vel);
    vel.rhs->add_vector(Fe, dof_indices_vel);
  }
}
