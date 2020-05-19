////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

void netfcfvfe::VelAssembly::assemble() {

  // Get required system alias
  // auto &vel = d_model_p->get_vel_assembly();
  auto &pres = d_model_p->get_pres_assembly();
  auto &tum = d_model_p->get_tum_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const double factor_p = deck.d_assembly_factor_p_t;

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_pres_neigh;

  // Store current and old solution
  Real chem_tum_cur = 0.;
  Real pres_cur = 0.;

  Gradient pres_grad;
  Gradient tum_grad;

  // Store current and old solution of neighboring element
  Real pres_neigh_cur = 0.;

  Gradient compute_rhs = 0.;

  auto side_n = std::vector<unsigned int>(deck.d_dim, 0);

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);
    pres.init_dof(elem);

    // init fe and element matrix and vector
    init_fe(elem);

    // get finite-volume quantities
    pres_cur = pres.get_current_sol(0);

    // Compute grad(p), which is constant in the element, at the center of
    // the element
    {
      pres_grad = 0.;
      for (auto &n : side_n)
        n = 0;

      // loop over sides of the element
      const auto elem_center = elem->centroid();
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);
          const auto neighbor_center = neighbor->centroid();
          const auto dx = elem_center - neighbor_center;

          // get dof id
          pres.init_dof(neighbor, dof_indices_pres_neigh);
          pres_neigh_cur = pres.get_current_sol(0, dof_indices_pres_neigh);

          // check if dx is aligned along x or y or z direction
          std::string dx_align = util::get_vec_major_axis(dx);

          if (dx_align == "x") {
            side_n[0]++;
            pres_grad(0) += (pres_cur - pres_neigh_cur) / dx(0);
          } else if (dx_align == "y") {
            side_n[1]++;
            pres_grad(1) += (pres_cur - pres_neigh_cur) / dx(1);
          } else if (dx_align == "z" and deck.d_dim > 2) {
            side_n[2]++;
            pres_grad(2) += (pres_cur - pres_neigh_cur) / dx(2);
          }
        }
      } // loop over neighbors

      // divide by number of neighbors
      for (unsigned int i = 0; i < deck.d_dim; i++)
        pres_grad(i) = pres_grad(i) / side_n[i];
    }

    // loop over quadrature points
    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // Computing solution
      chem_tum_cur = 0.;
      tum_grad = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {

        chem_tum_cur += d_phi[l][qp] * tum.get_current_sol_var(l, 1);
        tum_grad.add_scaled(d_dphi[l][qp], tum.get_current_sol_var(l, 0));
      }

      compute_rhs =
          d_JxW[qp] * factor_p *
          (-deck.d_tissue_flow_coeff * (pres_grad - chem_tum_cur * tum_grad));

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        for (unsigned int d=0; d<deck.d_dim; d++)
          d_Fe_var[d](i) += compute_rhs(d) * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          for (unsigned int d=0; d<deck.d_dim; d++)
            d_Ke_var[d][d](i, j) +=
                d_JxW[qp] * factor_p * d_phi[j][qp] * d_phi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}
