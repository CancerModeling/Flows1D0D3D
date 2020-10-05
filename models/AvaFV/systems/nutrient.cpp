////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_nut_source(const std::string &test_name, const Point &x,
                      const std::vector<double> &x0, const double &r) {

  if (test_name != "test_tum_2")
    return 0.;

  // For cylinder, take z-coord as 0
  const Point xc = util::to_point({x0[0], x0[1], 0.});

  // spherical source
  Point dx = Point(x(0), x(1), 0.) - xc;
  if (dx.norm() < r)
    return 1.;

  return 0.;
}
} // namespace

Number avafv::initial_condition_nut(const Point &p, const Parameters &es,
                                    const std::string &system_name,
                                    const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void avafv::boundary_condition_nut(EquationSystems &es) {

  const auto *deck = es.parameters.get<InpDeck *>("input_deck");

  std::set<boundary_id_type> ids;
  if (deck->d_nutrient_bc_north)
    ids.insert(2);
  if (deck->d_nutrient_bc_south)
    ids.insert(0);
  if (deck->d_nutrient_bc_east)
    ids.insert(1);
  if (deck->d_nutrient_bc_west)
    ids.insert(3);

  // get variable ids
  auto &sys = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("nutrient"));

  // create constant function for bc and apply
  ConstFunction<Number> one(1);
  DirichletBoundary diri_bc(ids, vars, &one);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

// Assembly class
void avafv::NutAssembly::assemble() {
  assemble_face();
  assemble_1();
}

void avafv::NutAssembly::assemble_face() {
  avafv::assemble_diffusion(*this, this->d_model_p);
}

void avafv::NutAssembly::assemble_1() {

  // Get required system alias
  //auto &nut = d_model_p->get_nut_assembly();
  auto &tum = d_model_p->get_tum_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nec = d_model_p->get_nec_assembly();
  auto &taf = d_model_p->get_taf_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // local matrix and vector
  DenseMatrix<Number> Ke(1, 1);
  DenseVector<Number> Fe(1);

  // Store current and old solution
  Real nut_old = 0.;
  Real tum_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;

  Real tum_proj = 0.;
  Real hyp_proj = 0.;
  Real nec_proj = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    tum.init_dof(elem);
    hyp.init_dof(elem);
    nec.init_dof(elem);
    taf.init_dof(elem);

    // reset matrix and force
    Ke(0, 0) = 0.;
    Fe(0) = 0.;

    // get solution in this element
    nut_old = get_old_sol(0);
    tum_cur = tum.get_current_sol_var(0, 0);
    hyp_cur = hyp.get_current_sol(0);
    nec_cur = nec.get_current_sol(0);

    if (deck.d_assembly_method == 1) {

      compute_rhs = deck.d_elem_size * nut_old;

      compute_mat =
        deck.d_elem_size *
        (1. + dt * (deck.d_lambda_P * (tum_cur - hyp_cur - nec_cur) +
                    deck.d_lambda_Ph * hyp_cur));

    } else {

      tum_proj = util::project_concentration(tum_cur);
      hyp_proj = util::project_concentration(hyp_cur);
      nec_proj = util::project_concentration(nec_cur);

      compute_rhs = deck.d_elem_size * nut_old;

      compute_mat =
        deck.d_elem_size *
        (1. + dt * (deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) +
                    deck.d_lambda_Ph * hyp_proj));
    }

    // add artificial source if asked
    Real artificial_source =
      get_nut_source(deck.d_test_name, elem->centroid(),
                     deck.d_nut_source_center, deck.d_nut_source_radius) -
      nut_old;
    if (artificial_source > 0.)
      compute_rhs += deck.d_elem_size * dt * artificial_source;

    // add
    Ke(0, 0) += factor_nut * compute_mat;
    Fe(0) += factor_nut * compute_rhs;

    // add to matrix
    d_sys.matrix->add_matrix(Ke, d_dof_indices_sys, d_dof_indices_sys);

    // add to vector
    d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
  }

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}