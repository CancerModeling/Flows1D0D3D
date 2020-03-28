////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number avafv::initial_condition_taf(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"TAF");

  return 0.;
}

// Assembly class
void avafv::TafAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  assemble_vol();
}

void avafv::TafAssembly::assemble_face() {

  // call diffusion calculation function
  avafv::assemble_diffusion(d_model_p->get_taf_assembly());
}
void avafv::TafAssembly::assemble_vol() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &taf = d_model_p->get_ecm_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real taf_old = 0.;
  Real hyp_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    hyp.init_dof(elem);

    // const unsigned int n_dofs = taf.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get fields at this element
      taf_old = get_old_sol(0);

      // get projected values of species
      hyp_proj = util::project_concentration(hyp.get_current_sol(0));

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += taf_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_proj;
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_proj;
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