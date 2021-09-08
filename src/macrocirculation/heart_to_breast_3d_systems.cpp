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
  auto &pts = d_model_p->d_perf_pts;
  auto &ball_r = d_model_p->d_perf_ball_radii;
  auto &out_pres = d_model_p->d_perf_pres;
  auto &out_elems = d_model_p->d_perf_elems_3D;

  for (size_t I = 0; I < pts.size(); I++) {
    auto &out_fn_I = out_fns[I];
    const auto &pI = out_pres[I];
    const auto &RI = out_fn_I->d_r;
    const auto &xI = pts[I];

    // loop over elements
    for (const auto &elem_id : out_elems[I]) {

      const auto &elem = d_mesh.elem_ptr(elem_id);

      // init dof map
      init_dof(elem);

      // init fe
      init_fe(elem);

      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
        for (unsigned int i = 0; i < d_phi.size(); i++) {
          d_Fe(i) += d_JxW[qp] * input.d_Lp_art_cap * pI * (*out_fn_I)(d_qpoints[qp]) *d_phi[i][qp];

          for (unsigned int j = 0; j < d_phi.size(); j++)
            d_Ke(i, j) += d_JxW[qp] * input.d_Lp_art_cap * (*out_fn_I)(d_qpoints[qp]) *d_phi[i][qp] * d_phi[j][qp];
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
  const double dt = d_model_p->d_dt;
  const double t = d_model_p->d_time;
  auto &p_tis_field = d_model_p->d_p_tis;
  auto &K_field = d_model_p->d_K_cap_field;
  auto &Lp_field = d_model_p->d_Lp_cap_tis_field;
  std::vector<unsigned int> K_dof_indices;
  std::vector<unsigned int> Lp_dof_indices;

  double p_tis_cur = 0.;
  double lhs = 0.;
  double rhs = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_tis_field.init_dof(elem);
    K_field.get_dof_map().dof_indices(elem, K_dof_indices);
    Lp_field.get_dof_map().dof_indices(elem, Lp_dof_indices);

    // init fe
    init_fe(elem);

    // get K at this element
    double elem_K = K_field.current_solution(K_dof_indices[0]);
    double elem_Lp = Lp_field.current_solution(Lp_dof_indices[0]);

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
          d_Ke(i, j) += d_JxW[qp] * elem_K * d_dphi[j][qp] * d_dphi[i][qp] + lhs * d_phi[j][qp] * d_phi[i][qp];
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
  const double dt = d_model_p->d_dt;
  const double t = d_model_p->d_time;
  auto &p_cap_field = d_model_p->d_p_cap;
  auto &K_field = d_model_p->d_K_tis_field;
  auto &Lp_field = d_model_p->d_Lp_cap_tis_field;
  std::vector<unsigned int> K_dof_indices;
  std::vector<unsigned int> Lp_dof_indices;

  double p_cap_cur = 0.;
  double lhs = 0.;
  double rhs = 0.;

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    p_cap_field.init_dof(elem);
    K_field.get_dof_map().dof_indices(elem, K_dof_indices);
    Lp_field.get_dof_map().dof_indices(elem, Lp_dof_indices);

    // init fe
    init_fe(elem);

    // get K at this element
    double elem_K = K_field.current_solution(K_dof_indices[0]);
    double elem_Lp = Lp_field.current_solution(Lp_dof_indices[0]);

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
          d_Ke(i, j) += d_JxW[qp] * elem_K * d_dphi[j][qp] * d_dphi[i][qp] + lhs * d_phi[j][qp] * d_phi[i][qp];
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


} // namespace macrocirculation