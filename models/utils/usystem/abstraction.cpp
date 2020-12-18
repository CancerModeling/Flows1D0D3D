////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "abstraction.hpp"

double util::BaseAssembly::compute_qoi(const std::string &type, unsigned int local_var_id) {

  if (type == "linf")
    return d_sys.solution->linfty_norm();

  if (type != "mass" and type != "l1" and type != "l2")
    libmesh_error_msg("Error: Invalid type = " + type +
                      " for qoi calculation");

  Real qoi = 0.;
  Real cur_sol = 0.;
  Real cur_sol_l = 0.;
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);

    init_fe(elem);

    if (d_sys_name == "Pressure" or d_sys_name == "Nutrient") {
      cur_sol = get_current_sol(0);
      if (type == "mass" and cur_sol < 0.)
        cur_sol = 0.;

      if (type == "mass")
        qoi += elem->volume() * cur_sol;
      else if (type == "l1")
        qoi += elem->volume() * std::abs(cur_sol);
      else if (type == "l2")
        qoi += elem->volume() * cur_sol * cur_sol;
    } else {
      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
        cur_sol = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {

          if (d_num_vars == 1)
            cur_sol_l = get_current_sol(l);
          else
            cur_sol_l = get_current_sol_var(l, local_var_id);
          if (type == "mass" and cur_sol_l < 0.)
            cur_sol_l = 0.;

          cur_sol = d_phi[l][qp] * cur_sol_l;
        }

        if (type == "mass")
          qoi += d_JxW[qp] * cur_sol;
        else if (type == "l1")
          qoi += d_JxW[qp] * std::abs(cur_sol);
        else if (type == "l2")
          qoi += d_JxW[qp] * cur_sol * cur_sol;
      } // loop over quadrature points
    }
  } // loop over elements

  // communicate qoi from other processors
  Real total_qoi = 0.;
  MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  if (type == "l2")
    total_qoi = std::sqrt(total_qoi);

  return total_qoi;
}

void util::BaseAssembly::project_physical_range(unsigned int local_var_id, double min, double max) {

  // loop over nodes and modify dofs
  for (const auto &node : d_mesh.local_node_ptr_range()) {

    const auto &dof = node->dof_number(d_sys.number(), d_var_id[local_var_id], 0);

    auto val = d_sys.current_solution(dof);
    if (val < min)
      d_sys.solution->add(dof, -val + min);
    else if (val > max)
      d_sys.solution->add(dof, -val + max);
    else
      continue;
  }

  d_sys.solution->close();
  d_sys.update();
}