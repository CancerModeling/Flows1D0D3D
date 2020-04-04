////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_taf_source(const std::string &test_name, const Point &x, const
double &L) {

  if (test_name != "test_taf_2")
    return 0.;

  double L_source_x = 0.6 * L;
  double L_source_y = 0.4 * L;
  double L_source_z = 0.;
  const Point xc = Point(L_source_x, L_source_y, L_source_z);
  // spherical source
  if (false) {
    if ((x - xc).norm() < 0.1 * xc.norm())
      return 1.;
  } else {
    if (x(0) > 0.9 * L_source_x and x(0) < 1.1 * L_source_x and
        x(1) > 0.9 * L_source_y and x(1) < 1.1 * L_source_y) {

      //            compute_rhs += JxW[qp] * dt * deck->d_lambda_TAF *
      //                           std::sin(2. * (2. * M_PI) * x(2) / L_source);
      const Point x_plane = Point(x(0), x(1), 0.);
      const Point xc_plane = Point(L_source_x, L_source_y, 0.);
      if ((x_plane - xc_plane).norm() < 0.1 * xc.norm())
        return 1.;
    }
  }

  return 0.;
}
}

Number netfv::initial_condition_taf(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"TAF");

  return 0.;
}

// Assembly class
void netfv::TafAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfv::TafAssembly::assemble_face() {

  // call diffusion-advection calculation function
  if (d_model_p->get_input_deck().d_advection_active)
    netfv::assemble_diffusion_advection(d_model_p->get_taf_assembly(),
                                       d_model_p->get_pres_assembly(),
                                       d_model_p->get_tum_assembly());
  else
    netfv::assemble_diffusion(d_model_p->get_taf_assembly());
}

void netfv::TafAssembly::assemble_1() {

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
  Real hyp_cur = 0.;

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
      hyp_cur = hyp.get_current_sol(0);
      taf_old = get_old_sol(0);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += taf_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_cur;
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_cur;

      // add artificial source if asked
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF *
               get_taf_source(deck.d_test_name, elem->centroid(),
                              deck.d_domain_params[1]);
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

void netfv::TafAssembly::assemble_2() {

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

      // add artificial source if asked
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF *
               get_taf_source(deck.d_test_name, elem->centroid(),
                              deck.d_domain_params[1]);
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

void netfv::TafAssembly::assemble_3() {

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
  Real taf_cur = 0.;
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
      taf_cur = get_current_sol(0);

      // get projected values of species
      hyp_proj = util::project_concentration(hyp.get_current_sol(0));

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += taf_old * deck.d_elem_size;

      // source contribution (in assemble_2 this term was implicit and here
      // we consider it explicitly)
      Real a_source = deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_proj;
      Fe(0) += -a_source * taf_cur;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF * hyp_proj;

      // add artificial source if asked
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_TAF *
               get_taf_source(deck.d_test_name, elem->centroid(),
                              deck.d_domain_params[1]);
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