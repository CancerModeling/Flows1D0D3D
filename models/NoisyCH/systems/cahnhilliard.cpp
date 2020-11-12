////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "cahnhilliard.hpp"

#include <random>

namespace noisych {

Number initial_condition_cahnhilliard_random(const Point &p,
                                             const Parameters &es,
                                             const std::string &system_name,
                                             const std::string &var_name) {
  libmesh_assert_equal_to(system_name, "CahnHilliard");

  static std::default_random_engine generator(2 + es.get<unsigned int>("rank"));
  const auto mean = 0.63;
  const auto interval = 0.1;
  std::uniform_real_distribution<double> distribution(mean - interval, mean + interval);

  if (var_name == "concentration") {
    return distribution(generator);
  }

  return 0;
}

Number initial_condition_cahnhilliard_circle(const Point &p,
                                             const Parameters &es,
                                             const std::string &system_name,
                                             const std::string &var_name) {
  // circle with middle point (1/2,1/2) and radius a
  const auto a = 0.25;
  const auto dist_squared = std::pow(p(0) - 0.5, 2) + std::pow(p(1) - 0.5, 2) + std::pow(p(2) - 0.5, 2);

  if (dist_squared < a * a)
    return std::exp(1. - 1. / (1. - dist_squared / a * a));

  return 0;
}

CahnHilliardConfig::CahnHilliardConfig()
    : cubic_root_num_eigenfunctions(0),
      seed(0),
      scale(0),
      length(1),
      lower_bound(0),
      upper_bound(1),
      dt(1e-6),
      C_psi(0),
      epsilon(0),
      mobility_constant(0) {}

CahnHilliardConfig CahnHilliardConfig::from_parameter_file(const std::string &filename) {
  CahnHilliardConfig config;

  GetPot input(filename);

  config.cubic_root_num_eigenfunctions = input("cubic_root_num_eigenfunctions", config.cubic_root_num_eigenfunctions);
  config.seed = input("seed", config.seed);
  config.scale = input("scale", config.scale);
  config.length = input("length", config.length);
  config.lower_bound = input("lower_bound", config.lower_bound);
  config.upper_bound = input("upper_bound", config.upper_bound);
  config.dt = input("dt", config.dt);
  config.C_psi = input("C_psi", config.C_psi);
  config.epsilon = input("epsilon", config.epsilon);
  config.mobility_constant = input("mobility_constant", config.mobility_constant);

  return config;
}

CahnHilliardAssembly::CahnHilliardAssembly(const std::string &system_name,
                                           MeshBase &mesh,
                                           TransientLinearImplicitSystem &sys,
                                           const CahnHilliardConfig &config)
    : util::BaseAssembly(system_name, mesh, sys, 2, {sys.variable_number("concentration"), sys.variable_number("potential")}),
      d_noise_assembly(
        std::pow(config.cubic_root_num_eigenfunctions, 3),
        config.seed,
        config.scale,
        config.length,
        config.lower_bound,
        config.upper_bound),
      d_dt(config.dt),
      d_C_psi(config.C_psi),
      d_epsilon(config.epsilon),
      d_mobility_constant(config.mobility_constant) {}

void CahnHilliardAssembly::calculate_new_stochastic_coefficients(double dt) {
  d_noise_assembly.calculate_new_stochastic_coefficients(dt);
}

// Assembly class
void CahnHilliardAssembly::assemble() {
  assemble_1();
  d_noise_assembly.assemble(*this);
}

void CahnHilliardAssembly::assemble_1() {
  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {
    init_dof(elem);
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
      Real c_cur = 0.;
      Real c_old = 0.;

      for (unsigned int l = 0; l < d_phi.size(); l++) {
        c_old += d_phi[l][qp] * get_old_sol_var(l, 0);
        c_cur += d_phi[l][qp] * get_current_sol_var(l, 0);
      }

      const Real mobility = d_mobility_constant * pow(util::proj(c_old) * util::proj(1. - c_old), 2);
      // const Real mobility = 1.;

      const Real compute_rhs_c = d_JxW[qp] * (d_dt * c_old * (1 - c_old) + c_old);
      // const Real compute_rhs_c = d_JxW[qp] * c_old;
      const Real compute_rhs_mu = d_JxW[qp] * c_old * d_C_psi * (4 * std::pow(c_old, 2) - 6 * c_old - 1);

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        // concentration
        d_Fe_var[0](i) += compute_rhs_c * d_phi[i][qp];

        // potential
        d_Fe_var[1](i) += compute_rhs_mu * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          // concentration
          d_Ke_var[0][0](i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];

          // concentration-potential
          d_Ke_var[0][1](i, j) += d_JxW[qp] * d_dt * mobility * d_dphi[j][qp] * d_dphi[i][qp];

          // potential
          d_Ke_var[1][1](i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];

          // potential-concentration
          d_Ke_var[1][0](i, j) += -d_JxW[qp] * 3.0 * d_C_psi * d_phi[j][qp] * d_phi[i][qp];
          d_Ke_var[1][0](i, j) += -d_JxW[qp] * pow(d_epsilon, 2) * d_dphi[j][qp] * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe, d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}

} // namespace noisych