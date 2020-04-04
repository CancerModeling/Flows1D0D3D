////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

// Assembly class
void netfv::GradTafAssembly::solve() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Get required system alias
  // auto &grad_taf = d_model_p->get_grad_taf_assembly();
  auto &taf = d_model_p->get_taf_assembly();

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_taf_neigh;

  // Store current solution
  Real taf_cur = 0.;
  Real taf_neigh_cur = 0.;
  Gradient grad_taf_cur;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_var_dof(elem);
    taf.init_dof(elem);

    // compute grad(phi_TAF) which is constant in the element
    taf_cur = taf.get_current_sol(0);
    for (unsigned int i=0; i<dim; i++)
      grad_taf_cur(i) = 0.;
    {
      unsigned int nx = 0, ny = 0, nz = 0;
      // loop over sides of the element
      auto elem_center = elem->centroid();
      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);
          auto neighbor_center = neighbor->centroid();
          const auto dx = elem_center - neighbor_center;

          // get dof id
          taf.init_dof(neighbor, dof_indices_taf_neigh);
          taf_neigh_cur = taf.get_current_sol(0, dof_indices_taf_neigh);

          // check if dx is aligned along x or y or z direction
          std::string dx_align = util::get_vec_major_axis(dx);

          if (dx_align == "x") {
            nx++;
            grad_taf_cur(0) += (taf_cur - taf_neigh_cur) / dx(0);
          } else if (dx_align == "y") {
            ny++;
            grad_taf_cur(1) += (taf_cur - taf_neigh_cur) / dx(1);
          } else if (dx_align == "z") {
            nz++;
            grad_taf_cur(2) += (taf_cur - taf_neigh_cur) / dx(2);
          }
        }
      } // loop over neighbors

      // divide by number of neighbors
      grad_taf_cur(0) = grad_taf_cur(0) / nx;
      grad_taf_cur(1) = grad_taf_cur(1) / ny;
      if (dim > 2)
        grad_taf_cur(2) = grad_taf_cur(2) / nz;
    }

    // set solution at this element dof
    for (unsigned int i = 0; i < dim; i++)
      d_sys.solution->set(get_var_global_dof_id(0, i), grad_taf_cur(i));
  }

  d_sys.solution->close();
  d_sys.update();
}
