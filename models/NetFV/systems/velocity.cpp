////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

// Assembly class
void netfv::VelAssembly::solve() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Get required system alias
  // auto &vel = d_model_p->get_vel_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_pres_neigh;
  std::vector<unsigned int> dof_indices_tum_neigh;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var_neigh(2);

  // Store current and old solution
  Real tum_cur = 0.;
  Real chem_tum_cur = 0.;
  Real pres_cur = 0.;

  // Store current and old solution of neighboring element
  Real tum_neigh_cur = 0.;
  Real chem_tum_neigh_cur = 0.;
  Real pres_neigh_cur = 0.;

  Gradient velocity;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_var_dof(elem);
    pres.init_dof(elem);
    tum.init_var_dof(elem);

    // get solution at this element
    tum_cur = tum.get_current_sol_var(0, 0);
    chem_tum_cur = tum.get_current_sol_var(0, 1);
    pres_cur = pres.get_current_sol(0);;
    velocity = 0.;
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
          // pres
          pres.init_dof(neighbor, dof_indices_pres_neigh);
          pres_neigh_cur = pres.get_current_sol(0, dof_indices_pres_neigh);

          // tum
          tum.init_var_dof(neighbor, dof_indices_tum_neigh,
                           dof_indices_tum_var_neigh);
          tum_neigh_cur = tum.get_current_sol_var(0, 0,
                                                  dof_indices_tum_var_neigh);
          chem_tum_neigh_cur = tum.get_current_sol_var(0, 1,
                                                       dof_indices_tum_var_neigh);

          // check if dx is aligned along x or y or z direction
          std::string dx_align = util::get_vec_major_axis(dx);

          Real a = -deck.d_tissue_flow_coeff *
                   ((pres_cur - pres_neigh_cur) -
                    chem_tum_neigh_cur * (tum_cur - tum_neigh_cur));

          if (dx_align == "x") {
            nx++;
            velocity(0) += a / dx(0);
          } else if (dx_align == "y") {
            ny++;
            velocity(1) += a / dx(1);
          } else if (dx_align == "z") {
            nz++;
            velocity(2) += a / dx(2);
          }
        }
      } // loop over neighbors

      // divide by number of neighbors
      velocity(0) = velocity(0) / nx;
      velocity(1) = velocity(1) / ny;
      if (dim > 2)
        velocity(2) = velocity(2) / nz;
    }

    // set solution at this element dof
    d_sys.solution->set(get_var_global_dof_id(0, 0),
                        velocity(0));
    d_sys.solution->set(get_var_global_dof_id(0, 1),
                        velocity(1));
    if (dim > 2)
      d_sys.solution->set(get_var_global_dof_id(0, 2),
                          velocity(2));
  }

  d_sys.solution->close();
  d_sys.update();
}
