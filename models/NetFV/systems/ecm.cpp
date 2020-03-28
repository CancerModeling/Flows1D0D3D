////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

Number netfv::initial_condition_ecm(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name){

  libmesh_assert_equal_to(system_name,"ECM");

  if (var_name == "ecm") {

    const auto *deck = es.get<netfv::InputDeck *>("input_deck");
    const auto &ic_data = deck->d_ecm_ic_data;

    if (ic_data.d_type.empty())
      return 0.;
    else if (ic_data.d_type == "spherical") {

      Point dx = p - Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
          ic_data.d_geom_params[2]);
      double r = ic_data.d_geom_params[3];

      if (dx.norm() < r - 1.0E-12)
        return ic_data.d_val * util::exp_decay_function(dx.norm() / r, 4.);
      else
        return 0.;
    } else if (ic_data.d_type == "elliptical") {

      Point xc = Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
                       ic_data.d_geom_params[2]);
      std::vector<double> r = {ic_data.d_geom_params[3], ic_data
                               .d_geom_params[4], ic_data.d_geom_params[5]};
      const Point dx = p - xc;

      // transform ellipse into ball of radius
      double ball_r = 0.;
      for (unsigned int i = 0; i < deck->d_dim; i++)
        ball_r = r[i] * r[i];
      ball_r = std::sqrt(ball_r);

      Point p_ball = util::ellipse_to_ball(p, xc, r,
                                                 deck->d_dim, ball_r);

      if (p_ball.norm() < ball_r - 1.0E-12) {

        return ic_data.d_val * util::exp_decay_function(p_ball.norm() / ball_r, 4.);
      } else
        return 0.;

    } else if (ic_data.d_type == "box") {

      Point x1 = Point(ic_data.d_geom_params[0],ic_data.d_geom_params[1],
                           ic_data.d_geom_params[2]);
      Point x2 = Point(ic_data.d_geom_params[3],ic_data.d_geom_params[4],
                       ic_data.d_geom_params[5]);

      if (util::is_inside_box(p, {x1, x2}))
        return ic_data.d_val;
      else
        return 0.;
    } else if (ic_data.d_type == "constant") {
      return ic_data.d_val;
    }
  }

  return 0.;
}

// Assembly class
void netfv::EcmAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_face();

  if (deck.d_assembly_method == 1)
    assemble_1();
  else if (deck.d_assembly_method == 2)
    assemble_2();
  else if (deck.d_assembly_method == 3)
    assemble_3();
}

void netfv::EcmAssembly::assemble_face() {

  // call advection calculation function
  netfv::assemble_advection(d_model_p->get_ecm_assembly(),
                             d_model_p->get_pres_assembly(),
                             d_model_p->get_tum_assembly());
}

void netfv::EcmAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &ecm = d_model_p->get_ecm_assembly();
  auto &nut = d_model_p->get_nut_assembly();  
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real ecm_old = 0.;
  Real ecm_cur = 0.;
  Real nut_cur = 0.;
  Real mde_cur = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);    
    mde.init_dof(elem);

    // const unsigned int n_dofs = ecm.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get remaining fields at this element
      ecm_cur = get_current_sol(0);
      ecm_old = get_old_sol(0);
      mde_cur = mde.get_current_sol(0);
      nut_cur = nut.get_current_sol(0);

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += ecm_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_ECM_D * mde_cur +
                       deck.d_lambda_ECM_P * nut_cur *
                           util::heaviside(ecm_cur - deck.d_bar_phi_ECM_P));
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_ECM_P * nut_cur *
               util::heaviside(ecm_cur - deck.d_bar_phi_ECM_P);
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

void netfv::EcmAssembly::assemble_2() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &ecm = d_model_p->get_ecm_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real ecm_old = 0.;
  Real ecm_proj = 0.;
  Real nut_proj = 0.;
  Real mde_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    mde.init_dof(elem);

    // const unsigned int n_dofs = ecm.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get remaining fields at this element
      ecm_old = get_old_sol(0);

      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      ecm_proj = util::project_concentration(get_current_sol(0));
      mde_proj = util::project_concentration(mde.get_current_sol(0));

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += ecm_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_ECM_D * mde_proj +
                       deck.d_lambda_ECM_P * nut_proj *
                       util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));
      Ke(0,0) += a_source;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_ECM_P * nut_proj *
               util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P);
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

void netfv::EcmAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();

  // Get required system alias
  // auto &ecm = d_model_p->get_ecm_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &mde = d_model_p->get_mde_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");

  // local matrix and vector
  DenseMatrix<Number> Ke(1,1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real ecm_old = 0.;
  Real ecm_cur = 0.;
  Real ecm_proj = 0.;
  Real nut_proj = 0.;
  Real mde_proj = 0.;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    nut.init_dof(elem);
    mde.init_dof(elem);

    // const unsigned int n_dofs = ecm.d_dof_indices_sys.size();

    // reset matrix and force
    Ke(0,0) = 0.;
    Fe(0) = 0.;

    // volume terms
    {
      // get remaining fields at this element
      ecm_old = get_old_sol(0);
      ecm_cur = get_current_sol(0);

      // get projected values of species
      nut_proj = util::project_concentration(nut.get_current_sol(0));
      ecm_proj = util::project_concentration(get_current_sol(0));
      mde_proj = util::project_concentration(mde.get_current_sol(0));

      // mass matrix
      Ke(0,0) += deck.d_elem_size;

      // previous time step term
      Fe(0) += ecm_old * deck.d_elem_size;

      // matrix contribution
      Real a_source = deck.d_elem_size * dt *
                      (deck.d_lambda_ECM_D * mde_proj +
                       deck.d_lambda_ECM_P * nut_proj *
                       util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));
      Fe(0) += -a_source * ecm_cur;

      // add source
      Fe(0) += deck.d_elem_size * dt * deck.d_lambda_ECM_P * nut_proj *
               util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P);
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