////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

// Assembly class
void netfv::GradTafAssembly::solve() {

  // Get required system alias
  // auto &grad_taf = d_model_p->get_grad_taf_assembly();
  auto &taf = d_model_p->get_taf_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_taf_neigh;

  // Store current solution
  Real taf_cur = 0.;
  Real taf_neigh_cur = 0.;
  Gradient taf_grad;

  Point elem_center = Point();
  Point neighbor_center = Point();

  auto ns = std::vector<unsigned int>(deck.d_dim, 0);

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    taf.init_dof(elem);

    elem_center = elem->centroid();

    // compute grad(phi_TAF) which is constant in the element
    taf_cur = taf.get_current_sol(0);
    taf_grad = 0.;
    for (auto &n : ns)
      n = 0;

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);
        neighbor_center = neighbor->centroid();
        const auto dx = elem_center - neighbor_center;

        // get dof id
        taf.init_dof(neighbor, dof_indices_taf_neigh);
        taf_neigh_cur = taf.get_current_sol(0, dof_indices_taf_neigh);

        // check if dx is aligned along x or y or z direction
        std::string dx_align = util::get_vec_major_axis(dx);

        if (dx_align == "x") {
          ns[0]++;
          taf_grad(0) += (taf_cur - taf_neigh_cur) / dx(0);
        } else if (dx_align == "y") {
          ns[1]++;
          taf_grad(1) += (taf_cur - taf_neigh_cur) / dx(1);
        } else if (dx_align == "z" and deck.d_dim > 2) {
          ns[2]++;
          taf_grad(2) += (taf_cur - taf_neigh_cur) / dx(2);
        }
      }
    } // loop over neighbors

    // set solution (after dividing by neighbor in x,y,z axis
    for (unsigned int i = 0; i < deck.d_dim; i++)
      d_sys.solution->set(get_var_global_dof_id(0, i), taf_grad(i) / ns[i]);
  }

  d_sys.solution->close();
  d_sys.update();
}