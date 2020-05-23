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
      hx = xmax  / nelx;
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
                        const std::vector<Number> &sol, double scale =
                         1. ) {

  // Looping through elements
  for (const auto &elem : sys.d_mesh.active_local_element_ptr_range()) {

    sys.init_dof(elem);
    sys.d_sys.solution->set(sys.get_global_dof_id(0), scale * sol[elem->id()]);
  }

  sys.d_sys.solution->close();
  sys.d_sys.update();
}

} // namespace util

#endif // UTIL_MODELUTIL_H
