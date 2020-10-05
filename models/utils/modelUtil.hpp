////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_MODELUTIL_H
#define UTIL_MODELUTIL_H

#include "umodel/model.hpp"
#include "usystem/abstraction.hpp"

namespace util {

/*!
 * @brief Create mesh
 */
inline void create_mesh(InpDeck &input, ReplicatedMesh &mesh) {
  double xmax = input.d_domain_params[1], ymax = input.d_domain_params[3],
         zmax = input.d_domain_params[5];
  unsigned int nelx, nely, nelz;
  double hx, hy, hz;
  if (input.d_dim == 2) {

    if (input.d_use_mesh_size_for_disc) {
      hx = input.d_mesh_size;
      hy = input.d_mesh_size;
      nelx = xmax / hx;
      nely = ymax / hy;
    } else {
      nelx = input.d_num_elems;
      nely = input.d_num_elems;
      hx = xmax / nelx;
      hy = ymax / nely;

      input.d_mesh_size = hx;
    }

    input.d_elem_face_size = hx;
    input.d_elem_size = hx * hx;
    input.d_face_by_h = 1.;

    // check if length of element in x and y direction are same
    if (std::abs(hx - hy) > 0.001 * hx) {
      libmesh_error_msg("Size of element needs to be same in both direction, "
                        "ie. element needs to be square\n"
                        "If domain is rectangle than specify number of "
                        "elements in x and y different so that element is "
                        "square\n");
    }

    // either read or create mesh
    if (input.d_restart)
      mesh.read(input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_square(mesh, nelx, nely, 0., xmax, 0., ymax,
                                          QUAD4);

  } else if (input.d_dim == 3) {

    if (input.d_use_mesh_size_for_disc) {
      hx = input.d_mesh_size;
      hy = input.d_mesh_size;
      hz = input.d_mesh_size;
      nelx = xmax / hx;
      nely = ymax / hy;
      nelz = zmax / hz;
    } else {
      nelx = input.d_num_elems;
      nely = input.d_num_elems;
      nelz = input.d_num_elems;
      hx = xmax / nelx;
      hy = ymax / nely;
      hz = zmax / nelz;

      input.d_mesh_size = hx;
    }

    input.d_elem_face_size = hx * hx;
    input.d_elem_size = hx * hx * hx;
    input.d_face_by_h = hx;

    // check if length of element in x and y direction are same
    if (std::abs(hx - hy) > 0.001 * hx or std::abs(hx - hz) > 0.001 * hx) {
      libmesh_error_msg("Size of element needs to be same in all three "
                        "direction, ie. element needs to be square\n"
                        "If domain is cuboid than specify number of "
                        "elements in x, y and z different so that element is "
                        "cube\n");
    }

    // either read or create mesh
    if (input.d_restart)
      mesh.read(input.d_mesh_restart_file);
    else
      MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, 0., xmax, 0.,
                                        ymax, 0., zmax, HEX8);
  }
}

inline double compute_diff_qoi(const std::string &type, BaseAssembly &a,
                               BaseAssembly &b, unsigned int local_var_id_a = 0,
                               unsigned int local_var_id_b = 0) {

  Real qoi = 0.;
  if (type == "inf")
    qoi = DBL_MAX;

  Real cur_sol_a = 0.;
  Real cur_sol_b = 0.;
  Real cur_sol_a_l = 0.;
  Real cur_sol_b_l = 0.;
  for (const auto &elem : a.d_mesh.active_local_element_ptr_range()) {

    a.init_dof(elem);
    b.init_dof(elem);
    for (unsigned int qp = 0; qp < a.d_qrule.n_points(); qp++) {
      cur_sol_a = 0.;
      cur_sol_b = 0.;
      for (unsigned int l = 0; l < a.d_phi.size(); l++) {

        if (a.d_num_vars == 1)
          cur_sol_a_l = a.get_current_sol(l);
        else
          cur_sol_a_l = a.get_current_sol_var(l, local_var_id_a);

        if (b.d_num_vars == 1)
          cur_sol_b_l = b.get_current_sol(l);
        else
          cur_sol_b_l = b.get_current_sol_var(l, local_var_id_b);

        if (type == "mass" and cur_sol_a_l < 0.)
          cur_sol_a_l = 0.;
        if (type == "mass" and cur_sol_b_l < 0.)
          cur_sol_b_l = 0.;

        cur_sol_a = a.d_phi[l][qp] * cur_sol_a_l;
        cur_sol_b = a.d_phi[l][qp] * cur_sol_b_l;
      }

      if (type == "mass")
        qoi += a.d_JxW[qp] * (cur_sol_a - cur_sol_b);
      else if (type == "l1")
        qoi += a.d_JxW[qp] * std::abs(cur_sol_a - cur_sol_b);
      else if (type == "l2")
        qoi += a.d_JxW[qp] * std::pow(cur_sol_a - cur_sol_b, 2);
      else if (type == "linf") {
        if (qoi > std::abs(cur_sol_a - cur_sol_b))
          qoi = std::abs(cur_sol_a - cur_sol_b);
      } else {
        libmesh_error_msg("Error: Invalid type = " + type + " for qoi calculation");
      }
    } // loop over quadrature points
  }   // loop over elements

  // communicate qoi from other processors
  Real total_qoi = 0.;
  if (type == "linf")
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  else
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

  if (type == "l2")
    total_qoi = std::sqrt(total_qoi);

  return total_qoi;
}

inline double compute_prolific_qoi(const std::string &type, BaseAssembly &a,
                                   BaseAssembly &b, BaseAssembly &c,
                                   unsigned int local_var_id_a = 0,
                                   unsigned int local_var_id_b = 0,
                                   unsigned int local_var_id_c = 0) {

  // a - tumor, b - hypoxic, c - necrotic

  Real qoi = 0.;
  if (type == "inf")
    qoi = DBL_MAX;

  Real cur_sol_a = 0.;
  Real cur_sol_b = 0.;
  Real cur_sol_c = 0.;
  Real cur_sol_a_l = 0.;
  Real cur_sol_b_l = 0.;
  Real cur_sol_c_l = 0.;
  for (const auto &elem : a.d_mesh.active_local_element_ptr_range()) {

    a.init_dof(elem);
    b.init_dof(elem);
    c.init_dof(elem);
    for (unsigned int qp = 0; qp < a.d_qrule.n_points(); qp++) {
      cur_sol_a = 0.;
      cur_sol_b = 0.;
      cur_sol_c = 0.;
      for (unsigned int l = 0; l < a.d_phi.size(); l++) {

        if (a.d_num_vars == 1)
          cur_sol_a_l = a.get_current_sol(l);
        else
          cur_sol_a_l = a.get_current_sol_var(l, local_var_id_a);

        if (b.d_num_vars == 1)
          cur_sol_b_l = b.get_current_sol(l);
        else
          cur_sol_b_l = b.get_current_sol_var(l, local_var_id_b);

        if (c.d_num_vars == 1)
          cur_sol_c_l = c.get_current_sol(l);
        else
          cur_sol_c_l = c.get_current_sol_var(l, local_var_id_c);

        if (type == "mass" and cur_sol_a_l < 0.)
          cur_sol_a_l = 0.;
        if (type == "mass" and cur_sol_b_l < 0.)
          cur_sol_b_l = 0.;
        if (type == "mass" and cur_sol_c_l < 0.)
          cur_sol_c_l = 0.;

        cur_sol_a = a.d_phi[l][qp] * cur_sol_a_l;
        cur_sol_b = a.d_phi[l][qp] * cur_sol_b_l;
        cur_sol_c = a.d_phi[l][qp] * cur_sol_c_l;
      }

      if (type == "mass")
        qoi += a.d_JxW[qp] * (cur_sol_a - cur_sol_b - cur_sol_c);
      else if (type == "l1")
        qoi += a.d_JxW[qp] * std::abs(cur_sol_a - cur_sol_b - cur_sol_c);
      else if (type == "l2")
        qoi += a.d_JxW[qp] * std::pow(cur_sol_a - cur_sol_b - cur_sol_c, 2);
      else if (type == "linf") {
        if (qoi > std::abs(cur_sol_a - cur_sol_b - cur_sol_c))
          qoi = std::abs(cur_sol_a - cur_sol_b - cur_sol_c);
      } else {
        libmesh_error_msg("Error: Invalid type = " + type + " for qoi calculation");
      }
    } // loop over quadrature points
  }   // loop over elements

  // communicate qoi from other processors
  Real total_qoi = 0.;
  if (type == "linf")
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  else
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

  if (type == "l2")
    total_qoi = std::sqrt(total_qoi);

  return total_qoi;
}

inline double compute_tumor_qoi(const std::string &type, BaseAssembly &a,
                                BaseAssembly &b, BaseAssembly &c,
                                unsigned int local_var_id_a = 0,
                                unsigned int local_var_id_b = 0,
                                unsigned int local_var_id_c = 0) {

  // a - prolific, b - hypoxic, c - necrotic

  Real qoi = 0.;
  if (type == "inf")
    qoi = DBL_MAX;

  Real cur_sol_a = 0.;
  Real cur_sol_b = 0.;
  Real cur_sol_c = 0.;
  Real cur_sol_a_l = 0.;
  Real cur_sol_b_l = 0.;
  Real cur_sol_c_l = 0.;
  for (const auto &elem : a.d_mesh.active_local_element_ptr_range()) {

    a.init_dof(elem);
    b.init_dof(elem);
    c.init_dof(elem);
    for (unsigned int qp = 0; qp < a.d_qrule.n_points(); qp++) {
      cur_sol_a = 0.;
      cur_sol_b = 0.;
      cur_sol_c = 0.;
      for (unsigned int l = 0; l < a.d_phi.size(); l++) {

        if (a.d_num_vars == 1)
          cur_sol_a_l = a.get_current_sol(l);
        else
          cur_sol_a_l = a.get_current_sol_var(l, local_var_id_a);

        if (b.d_num_vars == 1)
          cur_sol_b_l = b.get_current_sol(l);
        else
          cur_sol_b_l = b.get_current_sol_var(l, local_var_id_b);

        if (c.d_num_vars == 1)
          cur_sol_c_l = c.get_current_sol(l);
        else
          cur_sol_c_l = c.get_current_sol_var(l, local_var_id_c);

        if (type == "mass" and cur_sol_a_l < 0.)
          cur_sol_a_l = 0.;
        if (type == "mass" and cur_sol_b_l < 0.)
          cur_sol_b_l = 0.;
        if (type == "mass" and cur_sol_c_l < 0.)
          cur_sol_c_l = 0.;

        cur_sol_a = a.d_phi[l][qp] * cur_sol_a_l;
        cur_sol_b = a.d_phi[l][qp] * cur_sol_b_l;
        cur_sol_c = a.d_phi[l][qp] * cur_sol_c_l;
      }

      if (type == "mass")
        qoi += a.d_JxW[qp] * (cur_sol_a + cur_sol_b + cur_sol_c);
      else if (type == "l1")
        qoi += a.d_JxW[qp] * std::abs(cur_sol_a + cur_sol_b + cur_sol_c);
      else if (type == "l2")
        qoi += a.d_JxW[qp] * std::pow(cur_sol_a + cur_sol_b + cur_sol_c, 2);
      else if (type == "linf") {
        if (qoi > std::abs(cur_sol_a + cur_sol_b + cur_sol_c))
          qoi = std::abs(cur_sol_a + cur_sol_b + cur_sol_c);
      } else {
        libmesh_error_msg("Error: Invalid type = " + type + " for qoi calculation");
      }
    } // loop over quadrature points
  }   // loop over elements

  // communicate qoi from other processors
  Real total_qoi = 0.;
  if (type == "linf")
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  else
    MPI_Allreduce(&qoi, &total_qoi, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

  if (type == "l2")
    total_qoi = std::sqrt(total_qoi);

  return total_qoi;
}

inline void scale_pres(const ReplicatedMesh &mesh, util::BaseAssembly &pres,
                       const double &scale, std::vector<Number> &p_save,
                       std::vector<unsigned int> &p_dofs) {

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    pres.init_dof(elem);
    double pt = pres.get_current_sol(0);

    p_save.push_back(pt);
    p_dofs.push_back(pres.get_global_dof_id(0));

    pt = pt / scale;

    pres.d_sys.solution->set(pres.get_global_dof_id(0), pt);
  }

  pres.d_sys.solution->close();
  pres.d_sys.update();
}

inline void store_pres(const ReplicatedMesh &mesh, util::BaseAssembly &pres,
                       const double &scale, const std::vector<Number> &p_save,
                       const std::vector<unsigned int> &p_dofs) {

  for (unsigned int i = 0; i < p_dofs.size(); i++) {

    pres.d_sys.solution->set(p_dofs[i], p_save[i]);
  }
  pres.d_sys.solution->close();
  pres.d_sys.update();
}

inline void set_elem_sol(util::BaseAssembly &sys,
                         const std::vector<double> &sol, double scale = 1.) {

  // Looping through elements
  for (const auto &elem : sys.d_mesh.active_local_element_ptr_range()) {

    sys.init_dof(elem);
    sys.d_sys.solution->set(sys.get_global_dof_id(0), scale * sol[elem->id()]);
  }

  sys.d_sys.solution->close();
  sys.d_sys.update();
}

inline void get_elem_sol(util::BaseAssembly &sys, std::vector<double> &sol,
                         double scale = 1.) {

  if (sol.size() < sys.d_mesh.n_elem())
    sol = std::vector<double>(sys.d_mesh.n_elem(), 0.);

  // Looping through elements
  for (const auto &elem : sys.d_mesh.active_local_element_ptr_range()) {

    sys.init_dof(elem);
    sol[elem->id()] = scale * sys.get_current_sol(0);
  }
}

inline void localize_solution_with_elem_id_numbering_const_elem(
  util::BaseAssembly &sys, std::vector<double> &collect_sol,
  std::vector<double> &localize_sol, std::vector<unsigned int> var_ids = {0},
  bool resize_vec = true) {

  // gather solution in all processors
  // sys.d_sys.current_local_solution->localize(collect_sol);
  sys.d_sys.solution->localize(collect_sol);

  // check if we need to resize the vector
  auto num_vars = var_ids.size();
  if (localize_sol.size() != sys.d_mesh.n_elem() * num_vars) {
    if (resize_vec)
      localize_sol.resize(sys.d_mesh.n_elem() * num_vars);
    else
      libmesh_error_msg(
        "localize_sol size should match collect_sol size for system " +
        sys.d_sys_name);
  }

  // check if system has only 1 variable
  bool has_multiple_vars = sys.d_num_vars > 1;
  if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
    libmesh_error_msg("Invalid var ids to collect and localize solution");

  for (const auto &elem : sys.d_mesh.active_element_ptr_range()) {

    sys.init_dof(elem);
    int counter = 0;
    for (auto var : var_ids) {

      // compute value at center
      auto val = 0.;
      if (has_multiple_vars)
        val += collect_sol[sys.get_var_global_dof_id(0, var)];
      else
        val += collect_sol[sys.get_global_dof_id(0)];

      localize_sol[elem->id() * num_vars + counter] = val;
      counter++;
    }
  }
}

inline void localize_solution_with_elem_id_numbering_const_elem(
  util::BaseAssembly &sys,
  std::unique_ptr<NumericVector<Number>> &collect_sol,
  std::vector<double> &localize_sol, std::vector<unsigned int> var_ids = {0},
  bool resize_vec = true) {

  // gather solution in all processors
  // sys.d_sys.current_local_solution->localize(collect_sol);
  sys.d_sys.solution->localize(*collect_sol);

  // check if we need to resize the vector
  auto num_vars = var_ids.size();
  if (localize_sol.size() != sys.d_mesh.n_elem() * num_vars) {
    if (resize_vec)
      localize_sol.resize(sys.d_mesh.n_elem() * num_vars);
    else
      libmesh_error_msg(
        "localize_sol size should match collect_sol size for system " +
        sys.d_sys_name);
  }

  // check if system has only 1 variable
  bool has_multiple_vars = sys.d_num_vars > 1;
  if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
    libmesh_error_msg("Invalid var ids to collect and localize solution");

  for (const auto &elem : sys.d_mesh.active_element_ptr_range()) {

    sys.init_dof(elem);
    int counter = 0;
    for (auto var : var_ids) {

      // compute value at center
      auto val = 0.;
      if (has_multiple_vars)
        val += (*collect_sol)(sys.get_var_global_dof_id(0, var));
      else
        val += (*collect_sol)(sys.get_global_dof_id(0));

      localize_sol[elem->id() * num_vars + counter] = val;
      counter++;
    }
  }
}

inline void localize_solution_with_elem_id_numbering_non_const_elem(
  util::BaseAssembly &sys, std::vector<double> &collect_sol,
  std::vector<double> &localize_sol, std::vector<unsigned int> var_ids = {0},
  bool resize_vec = true) {

  // gather solution in all processors
  // sys.d_sys.current_local_solution->localize(collect_sol);
  sys.d_sys.solution->localize(collect_sol);

  // check if we need to resize the vector
  auto num_vars = var_ids.size();
  if (localize_sol.size() != sys.d_mesh.n_elem() * num_vars) {
    if (resize_vec)
      localize_sol.resize(sys.d_mesh.n_elem() * num_vars);
    else
      libmesh_error_msg(
        "localize_sol size should match collect_sol size for system " +
        sys.d_sys_name);
  }

  // check if system has only 1 variable
  bool has_multiple_vars = sys.d_num_vars > 1;
  if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
    libmesh_error_msg("Invalid var ids to collect and localize solution");

  for (const auto &elem : sys.d_mesh.active_element_ptr_range()) {

    sys.init_dof(elem);
    int counter = 0;
    for (auto var : var_ids) {

      // compute value at center
      auto val = 0.;
      for (unsigned int l = 0; l < sys.d_phi.size(); l++) {
        if (has_multiple_vars)
          val += collect_sol[sys.get_var_global_dof_id(l, var)];
        else
          val += collect_sol[sys.get_global_dof_id(l)];
      }
      val = val / (double(sys.d_phi.size()));

      localize_sol[elem->id() * num_vars + counter] = val;

      counter++;
    }
  }
}

inline void localize_solution_with_elem_id_numbering_non_const_elem(
  util::BaseAssembly &sys,
  std::unique_ptr<NumericVector<Number>> &collect_sol,
  std::vector<double> &localize_sol, std::vector<unsigned int> var_ids = {0},
  bool resize_vec = true) {

  // gather solution in all processors
  // sys.d_sys.current_local_solution->localize(collect_sol);
  sys.d_sys.solution->localize(*collect_sol);

  // check if we need to resize the vector
  auto num_vars = var_ids.size();
  if (localize_sol.size() != sys.d_mesh.n_elem() * num_vars) {
    if (resize_vec)
      localize_sol.resize(sys.d_mesh.n_elem() * num_vars);
    else
      libmesh_error_msg(
        "localize_sol size should match collect_sol size for system " +
        sys.d_sys_name);
  }

  // check if system has only 1 variable
  bool has_multiple_vars = sys.d_num_vars > 1;
  if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
    libmesh_error_msg("Invalid var ids to collect and localize solution");

  for (const auto &elem : sys.d_mesh.active_element_ptr_range()) {

    sys.init_dof(elem);
    int counter = 0;
    for (auto var : var_ids) {

      // compute value at center
      auto val = 0.;
      for (unsigned int l = 0; l < sys.d_phi.size(); l++) {
        if (has_multiple_vars)
          val += (*collect_sol)(sys.get_var_global_dof_id(l, var));
        else
          val += (*collect_sol)(sys.get_global_dof_id(l));
      }
      val = val / (double(sys.d_phi.size()));

      localize_sol[elem->id() * num_vars + counter] = val;

      counter++;
    }
  }
}

} // namespace util

#endif // UTIL_MODELUTIL_H
