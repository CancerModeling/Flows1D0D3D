////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"

void avafv::assemble_diffusion(BaseAssembly &sys) {

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
