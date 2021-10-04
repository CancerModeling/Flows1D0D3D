////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_3d_systems.hpp"

// makes the definition of HeartToBreast3DSolver precise
// HeartToBreast3DSolver is forward declared in include file heart_to_breast_3d_systems.hpp
#include "heart_to_breast_3d_solver.hpp"

namespace macrocirculation {

void CapillaryPressure::assemble_1d() {

  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;
  auto &out_fns = d_model_p->d_perf_fns;
  auto &out_pres = d_model_p->d_perf_pres;
  auto &out_pres_vein = d_model_p->d_perf_pres_vein;
  auto &out_elems = d_model_p->d_perf_elems_3D;

  for (size_t I = 0; I < out_pres.size(); I++) {
    auto &out_fn_I = out_fns[I];
    const auto &pI = out_pres[I];
    const auto &pvI = out_pres_vein[I];

    // loop over elements
    for (const auto &elem_id : out_elems[I]) {

      const auto &elem = d_mesh.elem_ptr(elem_id);

      // init dof map
      init_dof(elem);

      // init fe
      init_fe(elem);

      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
        for (unsigned int i = 0; i < d_phi.size(); i++) {
          d_Fe(i) += d_JxW[qp] * (input.d_Lp_art_cap * pI + input.d_Lp_vein_cap * pvI) *
                     (*out_fn_I)(d_qpoints[qp]) *d_phi[i][qp];

          for (unsigned int j = 0; j < d_phi.size(); j++)
            d_Ke(i, j) += d_JxW[qp] * (input.d_Lp_art_cap + input.d_Lp_vein_cap) *
                          (*out_fn_I)(d_qpoints[qp]) * d_phi[i][qp] * d_phi[j][qp];
        }
      } // loop over quad points

      d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                       d_dof_indices_sys);
      d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
      d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
    } // elem loop
  }   // outlet loop
}

// define assembly functions
void CapillaryPressure::assemble() {
  // call assemble term for artery-capillary coupling
  assemble_1d();

  // assemble from rest and then close the matrix and rhs vector
  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;
  auto &p_tis_field = d_model_p->d_p_tis;
  auto &N_bar_field = d_model_p->d_N_bar_cap_field;
  auto &N_bar_surf_field = d_model_p->d_N_bar_surf_cap_field;
  std::vector<unsigned int> elem_dof_indices;

  double p_tis_cur = 0.;
  double lhs = 0.;
  double rhs = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_tis_field.init_dof(elem);

    // init fe
    init_fe(elem);

    // get K_cap and Lp_cap_tis at this element
    N_bar_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Lp = input.d_Lp_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);

    N_bar_surf_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_K = input.d_K_cap * N_bar_surf_field.current_solution(elem_dof_indices[0]);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // get tissue pressure
      p_tis_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {
        p_tis_cur += d_phi[l][qp] * p_tis_field.get_current_sol(l);
      }

      lhs = d_JxW[qp] * input.d_rho_cap * elem_Lp;
      rhs = d_JxW[qp] * input.d_rho_cap * elem_Lp * p_tis_cur;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++)
          d_Ke(i, j) += d_JxW[qp] * input.d_rho_cap * elem_K * d_dphi[j][qp] * d_dphi[i][qp]
                        + lhs * d_phi[j][qp] * d_phi[i][qp];
      }
    } // loop over quad points

    // dirichlet bc
    if (false) {
      // The penalty value.
      const double penalty = 1.e10;

      for (auto s : elem->side_index_range())
        if (elem->neighbor_ptr(s) == nullptr) {
          init_face_fe(elem, s);

          for (unsigned int qp = 0; qp < d_qrule_face.n_points(); qp++) {
            // Matrix contribution
            for (unsigned int i = 0; i < d_phi_face.size(); i++) {
              d_Fe(i) += penalty * d_JxW_face[qp] * d_phi_face[i][qp];
              for (unsigned int j = 0; j < d_phi_face.size(); j++)
                d_Ke(i, j) += penalty * d_JxW_face[qp] * d_phi_face[i][qp] * d_phi_face[j][qp];
            }
          }
        }
    }

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void TissuePressure::assemble() {
  // assemble and close the matrix and rhs vector
  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;
  auto &p_cap_field = d_model_p->d_p_cap;
  auto &N_bar_field = d_model_p->d_N_bar_cap_field;
  auto &N_bar_surf_field = d_model_p->d_N_bar_surf_cap_field;
  std::vector<unsigned int> elem_dof_indices;

  double p_cap_cur = 0.;
  double lhs = 0.;
  double rhs = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_cap_field.init_dof(elem);

    // init fe
    init_fe(elem);

    // get K_cap and Lp_cap_tis at this element
    N_bar_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Lp = input.d_Lp_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);

    N_bar_surf_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_K = input.d_K_cap * N_bar_surf_field.current_solution(elem_dof_indices[0]);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // get tissue pressure
      p_cap_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {
        p_cap_cur += d_phi[l][qp] * p_cap_field.get_current_sol(l);
      }

      lhs = d_JxW[qp] * input.d_rho_cap * elem_Lp;
      rhs = d_JxW[qp] * input.d_rho_cap * elem_Lp * p_cap_cur;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++)
          d_Ke(i, j) += d_JxW[qp] * input.d_rho_cap * elem_K * d_dphi[j][qp] * d_dphi[i][qp]
                        + lhs * d_phi[j][qp] * d_phi[i][qp];
      }
    } // loop over quad points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void CapillaryNutrient::assemble_1d() {

  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;
  auto &out_fns = d_model_p->d_perf_fns;
  auto &out_elems = d_model_p->d_perf_elems_3D;
  auto &out_pres = d_model_p->d_perf_pres;
  auto &out_pres_vein = d_model_p->d_perf_pres_vein;
  auto &out_nut = d_model_p->d_perf_nut;
  auto &out_nut_vein = d_model_p->d_perf_nut_vein;
  auto &p_cap_field = d_model_p->d_p_cap;

  const double dt = d_model_p->d_dt;

  for (size_t I = 0; I < out_pres.size(); I++) {
    auto &out_fn_I = out_fns[I];
    const auto &pI = out_pres[I];
    const auto &pvI = out_pres_vein[I];
    const auto &nutI = out_nut[I];
    const auto &nutvI = out_nut_vein[I];
    const auto &RI = out_fn_I->d_r;

    // loop over elements
    for (const auto &elem_id : out_elems[I]) {

      const auto &elem = d_mesh.elem_ptr(elem_id);

      // init dof map
      init_dof(elem);
      p_cap_field.init_dof(elem);

      // init fe
      init_fe(elem);

      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
        double p_cap_qp = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {
          p_cap_qp += d_phi[l][qp] * p_cap_field.get_current_sol(l);
        }

        // artery-capillary coupling
        if (pI - p_cap_qp > 0.) {
          // to rhs
          for (unsigned int i = 0; i < d_phi.size(); i++) {
            d_Fe(i) += d_JxW[qp] * dt * (1. - input.d_rnut_art_cap) * input.d_Lp_art_cap
                       * (pI - p_cap_qp) * (*out_fn_I)(d_qpoints[qp]) * nutI * d_phi[i][qp];
          }
        }
        else {
          // to lhs
          for (unsigned int i = 0; i < d_phi.size(); i++) {
            for (unsigned int j = 0; j < d_phi.size(); j++) {
              d_Ke(i, j) += -d_JxW[qp] * dt * (1. - input.d_rnut_art_cap) * input.d_Lp_art_cap
                * (pI - p_cap_qp) * (*out_fn_I)(d_qpoints[qp]) * d_phi[j][qp] *d_phi[i][qp];
            }
          }
        } // artery-capillary

        // vein-capillary coupling
        if (pvI - p_cap_qp > 0.) {
          // to rhs
          for (unsigned int i = 0; i < d_phi.size(); i++) {
            d_Fe(i) += d_JxW[qp] * dt * (1. - input.d_rnut_vein_cap) * input.d_Lp_vein_cap
              * (pvI - p_cap_qp) * (*out_fn_I)(d_qpoints[qp]) * nutvI * d_phi[i][qp];
          }
        } else {
          // to lhs
          for (unsigned int i = 0; i < d_phi.size(); i++) {
            for (unsigned int j = 0; j < d_phi.size(); j++) {
              d_Ke(i, j) += -d_JxW[qp] * dt * (1. - input.d_rnut_vein_cap) * input.d_Lp_vein_cap
                * (pvI - p_cap_qp) * (*out_fn_I)(d_qpoints[qp]) * d_phi[j][qp] *d_phi[i][qp];
            }
          }
        } // vein-capillary

      } // loop over quad points

      d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                       d_dof_indices_sys);
      d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
      d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
    } // elem loop
  }   // outlet loop
}

// define assembly functions
void CapillaryNutrient::assemble() {
  // call assemble term for artery-capillary coupling
  assemble_1d();

  // assemble from rest and then close the matrix and rhs vector
  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;

  auto &p_cap_field = d_model_p->d_p_cap;
  auto &p_tis_field = d_model_p->d_p_tis;
  auto &nut_tis_field = d_model_p->d_nut_tis;
  auto &N_bar_field = d_model_p->d_N_bar_cap_field;
  auto &N_bar_surf_field = d_model_p->d_N_bar_surf_cap_field;
  std::vector<unsigned int> elem_dof_indices;

  const double dt = d_model_p->d_dt;

  double p_tis_cur = 0.;
  double p_cap_cur = 0.;
  double p_diff = 0.;
  double nut_tis_cur = 0.;
  double nut_cap_old = 0.;
  lm::Gradient p_cap_grad_cur = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_tis_field.init_dof(elem);
    p_cap_field.init_dof(elem);
    nut_tis_field.init_dof(elem);

    // init fe
    init_fe(elem);

    // get parameters at this element
    N_bar_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Lp = input.d_Lp_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);
    double elem_Lnut = input.d_Lnut_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);

    N_bar_surf_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Dnut = input.d_Dnut_cap * N_bar_surf_field.current_solution(elem_dof_indices[0]);
    double elem_K = input.d_K_cap * N_bar_surf_field.current_solution(elem_dof_indices[0]);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // get current values
      p_tis_cur = 0.;
      p_cap_cur = 0.;
      nut_tis_cur = 0.;
      nut_cap_old = 0.;
      p_cap_grad_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {
        p_tis_cur += d_phi[l][qp] * p_tis_field.get_current_sol(l);
        p_cap_cur += d_phi[l][qp] * p_cap_field.get_current_sol(l);
        nut_tis_cur += d_phi[l][qp] * nut_tis_field.get_current_sol(l);
        nut_cap_old += d_phi[l][qp] * get_old_sol(l);
        p_cap_grad_cur.add_scaled(d_dphi[l][qp],
                                  p_cap_field.get_current_sol(l));
      }
      p_diff = p_tis_cur - p_cap_cur;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += d_JxW[qp] * (nut_cap_old + dt * elem_Lnut * nut_tis_cur) * d_phi[i][qp];
        if (p_diff > 0.)
          d_Fe(i) += d_JxW[qp] * dt * (1. - input.d_rnut_cap) * elem_Lp * p_diff * nut_tis_cur * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {
          d_Ke(i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];
          d_Ke(i, j) += d_JxW[qp] * dt * elem_Dnut * d_dphi[j][qp] * d_dphi[i][qp];
          d_Ke(i, j) += d_JxW[qp] * dt * elem_Lnut * d_phi[j][qp] * d_phi[i][qp];

          if (!(p_diff > 0.))
            d_Ke(i, j) += -d_JxW[qp] * dt * (1. - input.d_rnut_cap) * elem_Lp * p_diff * d_phi[j][qp] * d_phi[i][qp];

          d_Ke(i, j) += d_JxW[qp] * dt * d_phi[j][qp] * elem_K * p_cap_grad_cur * d_dphi[i][qp];
        }

      }
    } // loop over quad points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void TissueNutrient::assemble() {
  // assemble and close the matrix and rhs vector
  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;

  auto &p_cap_field = d_model_p->d_p_cap;
  auto &p_tis_field = d_model_p->d_p_tis;
  auto &nut_cap_field = d_model_p->d_nut_cap;
  auto &K_tis_field = d_model_p->d_K_tis_field;
  auto &Dnut_tis_field = d_model_p->d_Dnut_tis_field;
  auto &N_bar_field = d_model_p->d_N_bar_cap_field;
  std::vector<unsigned int> elem_dof_indices;

  const double dt = d_model_p->d_dt;

  double p_tis_cur = 0.;
  double p_cap_cur = 0.;
  double p_diff = 0.;
  double nut_cap_cur = 0.;
  double nut_tis_old = 0.;
  lm::Gradient p_tis_grad_cur = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_tis_field.init_dof(elem);
    p_cap_field.init_dof(elem);
    nut_cap_field.init_dof(elem);

    // init fe
    init_fe(elem);

    // get parameters at this element
    N_bar_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Lp = input.d_Lp_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);
    double elem_Lnut = input.d_Lnut_cap_tis * N_bar_field.current_solution(elem_dof_indices[0]);

    K_tis_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_K = K_tis_field.current_solution(elem_dof_indices[0]);

    Dnut_tis_field.get_dof_map().dof_indices(elem, elem_dof_indices);
    double elem_Dnut = Dnut_tis_field.current_solution(elem_dof_indices[0]);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // get current values
      p_tis_cur = 0.;
      p_cap_cur = 0.;
      nut_cap_cur = 0.;
      nut_tis_old = 0.;
      p_tis_grad_cur = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {
        p_tis_cur += d_phi[l][qp] * p_tis_field.get_current_sol(l);
        p_cap_cur += d_phi[l][qp] * p_cap_field.get_current_sol(l);
        nut_cap_cur += d_phi[l][qp] * nut_cap_field.get_current_sol(l);
        nut_tis_old += d_phi[l][qp] * get_old_sol(l);
        p_tis_grad_cur.add_scaled(d_dphi[l][qp],
                                  p_tis_field.get_current_sol(l));
      }
      p_diff = p_tis_cur - p_cap_cur;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += d_JxW[qp] * (nut_tis_old + dt * elem_Lnut * nut_cap_cur) * d_phi[i][qp];
        if (p_diff < 0.)
          d_Fe(i) += -d_JxW[qp] * dt * (1. - input.d_rnut_cap) * elem_Lp * p_diff * nut_cap_cur * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {
          d_Ke(i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];
          d_Ke(i, j) += d_JxW[qp] * dt * elem_Dnut * d_dphi[j][qp] * d_dphi[i][qp];
          d_Ke(i, j) += d_JxW[qp] * dt * elem_Lnut * d_phi[j][qp] * d_phi[i][qp];

          if (!(p_diff < 0.))
            d_Ke(i, j) += d_JxW[qp] * dt * (1. - input.d_rnut_cap) * elem_Lp * p_diff * d_phi[j][qp] * d_phi[i][qp];

          d_Ke(i, j) += d_JxW[qp] * dt * d_phi[j][qp] * elem_K * p_tis_grad_cur * d_dphi[i][qp];
        }

      }
    } // loop over quad points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void Tumor::assemble() {
  // assemble and close the matrix and rhs vector
  auto &eq_sys = d_model_p->get_system();
  const auto &input = d_model_p->d_input;

  auto &nut_tis_field = d_model_p->d_nut_tis;
  const double dt = d_model_p->d_dt;

  double nut_tis_cur = 0.;
  double tum_old = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    nut_tis_field.init_dof(elem);

    // init fe
    init_fe(elem);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

      // get current values
      nut_tis_cur = 0.;
      tum_old = 0.;
      for (unsigned int l = 0; l < d_phi.size(); l++) {
        nut_tis_cur += d_phi[l][qp] * nut_tis_field.get_current_sol(l);
        tum_old += d_phi[l][qp] * get_old_sol_var(l, 0);
      }

      double mobility = input.d_tum_mob * std::pow(tum_old*(1. - tum_old), 2) + 1.e-8;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {

        d_Fe_var[0](i) += d_JxW[qp] * tum_old * d_phi[i][qp];

        d_Fe_var[1](i) += d_JxW[qp] * input.d_tum_dw * (4. * std::pow(tum_old, 3) - 6 * std::pow(tum_old, 2) - tum_old) * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++) {

          d_Ke_var[0][0](i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];
          d_Ke_var[0][0](i, j) += d_JxW[qp] * dt * (input.d_lambda_A - input.d_lambda_P * (1. - tum_old) * nut_tis_cur) * d_phi[j][qp] * d_phi[i][qp];

          d_Ke_var[0][1](i, j) += d_JxW[qp] * dt * mobility *  d_dphi[j][qp] * d_dphi[i][qp];

          d_Ke_var[1][1](i, j) += d_JxW[qp] * d_phi[j][qp] * d_phi[i][qp];

          d_Ke_var[1][0](i, j) += d_JxW[qp] * (-input.d_tum_dw * 3.) * d_phi[j][qp] * d_phi[i][qp];

          d_Ke_var[1][0](i, j) += d_JxW[qp] * (-input.d_tum_eps * input.d_tum_eps) * d_dphi[j][qp] * d_dphi[i][qp];
        }

      }
    } // loop over quad points

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

} // namespace macrocirculation