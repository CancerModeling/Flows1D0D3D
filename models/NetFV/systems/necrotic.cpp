////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfv::initial_condition_nec(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"Necrotic");

  return 0.;
}

// Assembly class
void netfv::NecAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2 || deck.d_assembly_method == 3)
    assemble_2();
}

void netfv::NecAssembly::assemble_face() {

  // call advection calculation function
  if (d_model_p->get_input_deck().d_advection_active)
    netfv::assemble_advection(d_model_p->get_nec_assembly(),
                              d_model_p->get_pres_assembly(),
                              d_model_p->get_tum_assembly());
}

void netfv::NecAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &nec = d_model_p->get_nec_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  
  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nec_old = 0.;
  Real nut_cur = 0.;
  Real hyp_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    hyp.init_dof(elem);
    
    // const unsigned int n_dofs = nec.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get solution in this element
      nut_cur = nut.get_current_sol(0);
      nec_old = get_old_sol(0);
      hyp_cur = hyp.get_current_sol(0);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += nec_old * deck.d_elem_size;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_HN * hyp_cur * util::heaviside
          (deck.d_sigma_HN - nut_cur);
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

void netfv::NecAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &nec = d_model_p->get_nec_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nec_old = 0.;
  Real nut_proj = 0.;
  Real hyp_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    hyp.init_dof(elem);

    // const unsigned int n_dofs = nec.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get solution in this element
      nec_old = get_old_sol(0);

      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      hyp_proj = util::project_concentration(hyp.get_current_sol(0));

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += nec_old * deck.d_elem_size;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_HN * hyp_proj * util::heaviside
          (deck.d_sigma_HN - nut_proj);
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