////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "cahnhilliard.hpp"
#include "../model.hpp"

#include <random>

namespace noisych {


Number initial_condition_cahnhilliard(const Point &p,
                                      const Parameters &es,
                                      const std::string &system_name,
                                      const std::string &var_name) {
  libmesh_assert_equal_to(system_name, "CahnHilliard");

  static std::default_random_engine generator( 2 + es.get<unsigned int>("rank") );
  const auto mean = 0.63;
  const auto interval = 0.1;
  std::uniform_real_distribution<double> distribution(mean-interval,mean+interval);

  if (var_name == "concentration") {
    return distribution(generator);
  }

  return 0;
}

CahnHilliardAssembly::CahnHilliardAssembly(const std::string &system_name,
                                           MeshBase &mesh,
                                           TransientLinearImplicitSystem &sys)
    : util::BaseAssembly(system_name, mesh, sys, 2, {sys.variable_number("concentration"), sys.variable_number("potential")}),
      d_noise_assembly(
        5*5*5,
        42,
        0.05,
        1.,
        0.4,
        0.6),
      d_dt(5e-6),
      d_C_psi(100),
      d_epsilon(0.1),
      d_mobility_constant(32.)
{}

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
          d_Ke_var[1][0](i, j) += - d_JxW[qp] * 3.0 * d_C_psi * d_phi[j][qp] * d_phi[i][qp];
          d_Ke_var[1][0](i, j) += - d_JxW[qp] * pow(d_epsilon, 2) * d_dphi[j][qp] * d_dphi[i][qp];
        }
      }
    } // loop over quadrature points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe, d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  }
}

} // namespace noisych