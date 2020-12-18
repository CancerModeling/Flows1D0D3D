////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_SYSTEM_ABSTRACTION_H
#define UTIL_SYSTEM_ABSTRACTION_H

#include "utils.hpp"

namespace util {

/*!
 * @brief Class to perform assembly of tumor species
 */
class BaseAssembly : public System::Assembly {

public:
  /*!
   * @brief Constructor
   *
   * @param sys_name Name of system
   * @param sys System
   */
  BaseAssembly(const std::string system_name, MeshBase &mesh,
               TransientLinearImplicitSystem &sys, const unsigned int &num_vars,
               const std::vector<unsigned int> &var_id)
      : d_sys_name(system_name), d_mesh(mesh), d_sys(sys), d_num_vars(num_vars),
        d_var_id(var_id), d_n_dofs(0), d_n_dofs_var(0),
        d_dof_map_sys(d_sys.get_dof_map()), d_fe_type(d_sys.variable_type(0)),
        d_qrule(d_mesh.mesh_dimension(), d_fe_type.default_quadrature_order()),
        d_fe(FEBase::build(d_mesh.mesh_dimension(), d_fe_type)),
        d_qpoints(d_fe->get_xyz()), d_JxW(d_fe->get_JxW()),
        d_phi(d_fe->get_phi()), d_dphi(d_fe->get_dphi()),
        d_qrule_face(d_mesh.mesh_dimension() - 1,
                     d_fe_type.default_quadrature_order()),
        d_fe_face(FEBase::build(d_mesh.mesh_dimension(), d_fe_type)),
        d_qpoints_face(d_fe_face->get_xyz()), d_JxW_face(d_fe_face->get_JxW()),
        d_phi_face(d_fe_face->get_phi()), d_dphi_face(d_fe_face->get_dphi()),
        d_qface_normals(d_fe_face->get_normals()), d_localized_sol(nullptr),
        d_assemble_matrix(true), d_implicit_assembly(true) {

    d_dof_indices_sys_var.resize(d_num_vars);

    d_fe->attach_quadrature_rule(&d_qrule);
    d_fe_face->attach_quadrature_rule(&d_qrule_face);

    if (d_num_vars > 1) {

      d_Fe_var = std::vector<DenseSubVector<Number>>(d_num_vars, d_Fe);
      d_Ke_var = std::vector<std::vector<DenseSubMatrix<Number>>>(
        d_num_vars, std::vector<DenseSubMatrix<Number>>(d_num_vars, d_Ke));
    }
  }

  /*!
   * @brief Assembly function
   *
   * To be overridden by inheriting classes
   */
  virtual void assemble() = 0;

  /*!
   * @brief Calls system solver
   *
   * This function is virtual for flexibility. Default method should work.
   */
  virtual void solve() { d_sys.solve(); };

  /*!
   * @brief Calls custom solver
   *
   * This function is virtual for flexibility. Default method should work.
   */
  virtual void solve_custom(){};

  /*!
   * @brief compute norm of the system
   */
  double compute_qoi(const std::string &type, unsigned int local_var_id = 0);

  /*!
   * @brief project solution to the physical range
   */
  void project_physical_range(unsigned int local_var_id = 0, double min = 0., double max = 1.);

  /*!
   * @brief Initializes the fe and local matrix and vector
   *
   * @param elem Pointer to the element
   */
  void init_fe(const Elem *elem) {

    d_fe->reinit(elem);

    d_Fe.resize(d_n_dofs);
    d_Ke.resize(d_n_dofs, d_n_dofs);

    if (d_num_vars > 1)
      for (unsigned int i = 0; i < d_num_vars; i++) {

        d_Fe_var[i].reposition(i * d_n_dofs_var, d_n_dofs_var);

        for (unsigned int j = 0; j < d_num_vars; j++)
          d_Ke_var[i][j].reposition(i * d_n_dofs_var, j * d_n_dofs_var,
                                    d_n_dofs_var, d_n_dofs_var);
      }
  }

  /*!
   * @brief Initializes the fe and local matrix and vector
   *
   * @param elem Pointer to the element
   */
  void init_face_fe(const Elem *elem, unsigned short side) {
    d_fe_face->reinit(elem, side);
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * @param elem Pointer to the element
   */
  void init_dof(const Elem *elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
    d_n_dofs = d_dof_indices_sys.size();
    d_n_dofs_var = d_n_dofs;

    if (d_num_vars > 1) {
      for (unsigned int var = 0; var < d_num_vars; var++)
        d_dof_map_sys.dof_indices(elem, d_dof_indices_sys_var[var],
                                  d_var_id[var]);

      d_n_dofs_var = d_dof_indices_sys_var[0].size();
    }
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * For system with more than one variable
   *
   * @param elem Pointer to the element
   */
  void init_var_dof(const Elem *elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
    for (unsigned int var = 0; var < d_num_vars; var++)
      d_dof_map_sys.dof_indices(elem, d_dof_indices_sys_var[var],
                                d_var_id[var]);

    d_n_dofs = d_dof_indices_sys.size();
    d_n_dofs_var = d_dof_indices_sys_var[0].size();
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * It populates the dof_indices passed to the function
   *
   * @param elem Pointer to the element
   * @param dof_indices_sys Dof indices vector
   */
  void init_dof(const Elem *elem, std::vector<unsigned int> &dof_indices_sys) {

    d_dof_map_sys.dof_indices(elem, dof_indices_sys);
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * For system with more than one variable
   *
   * It populates the dof_indices passed to the function
   *
   * @param elem Pointer to the element
   * @param dof_indices_sys Dof indices vector
   * @param dof_indices_sys_var Vector of vector of dof indices
   */
  void
  init_var_dof(const Elem *elem, std::vector<unsigned int> &dof_indices_sys,
               std::vector<std::vector<unsigned int>> &dof_indices_sys_var) {

    d_dof_map_sys.dof_indices(elem, dof_indices_sys);
    for (unsigned int var = 0; var < d_num_vars; var++)
      d_dof_map_sys.dof_indices(elem, dof_indices_sys_var[var], d_var_id[var]);
  }

  /*!
   * @brief Returns global dof id
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @return id Global dof id
   */
  unsigned int get_global_dof_id(const unsigned int &local_dof_id) {
    return d_dof_indices_sys[local_dof_id];
  }

  /*!
   * @brief Returns global dof id
   *
   * For system with more than one variable
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param local_var_id Local variable id (0 for tumor, 1 for chemical
   * potential)
   * @return id Global dof id
   */
  unsigned int get_var_global_dof_id(const unsigned int &local_dof_id,
                                     const unsigned int &local_var_id) {
    return d_dof_indices_sys_var[local_var_id][local_dof_id];
  }

  /*!
   * @brief Reads current solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @return solution Solution at the given dof
   */
  double get_current_sol(const unsigned int &local_dof_id) {
    return d_sys.current_solution(d_dof_indices_sys[local_dof_id]);
  }

  /*!
   * @brief Reads current solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param local_var_id Local variable id (0 for tumor, 1 for chemical
   * potential)
   * @return solution Solution at the given dof
   */
  double get_current_sol_var(const unsigned int &local_dof_id,
                             const unsigned int &local_var_id) {
    return d_sys.current_solution(
      d_dof_indices_sys_var[local_var_id][local_dof_id]);
  }

  /*!
   * @brief Reads current solution
   *
   * This function uses the supplied dof indices vector to compute solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param dof_indices_sys Dof indices vector
   * @return solution Solution at the given dof
   */
  double get_current_sol(const unsigned int &local_dof_id,
                         const std::vector<unsigned int> &dof_indices_sys) {
    return d_sys.current_solution(dof_indices_sys[local_dof_id]);
  }

  /*!
   * @brief Reads current solution
   *
   * This function uses the supplied dof indices vector to compute solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param local_var_id Local variable id (0 for tumor, 1 for chemical
   * potential)
   * @param dof_indices_sys_var Vector of Dof indices vector
   * @return solution Solution at the given dof
   */
  double get_current_sol_var(
    const unsigned int &local_dof_id, const unsigned int &local_var_id,
    const std::vector<std::vector<unsigned int>> &dof_indices_sys_var) {
    return d_sys.current_solution(
      dof_indices_sys_var[local_var_id][local_dof_id]);
  }

  /*!
   * @brief Reads old solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @return solution Solution at the given dof
   */
  double get_old_sol(const unsigned int &local_dof_id) {
    return d_sys.old_solution(d_dof_indices_sys[local_dof_id]);
  }

  /*!
   * @brief Reads old solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param local_var_id Local variable id (0 for tumor, 1 for chemical
   * potential)
   * @return solution Solution at the given dof
   */
  double get_old_sol_var(const unsigned int &local_dof_id,
                         const unsigned int &local_var_id) {
    return d_sys.old_solution(
      d_dof_indices_sys_var[local_var_id][local_dof_id]);
  }

  /*!
   * @brief Reads old solution
   *
   * This function uses the supplied dof indices vector to compute solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param dof_indices_sys Dof indices vector
   * @return solution Solution at the given dof
   */
  double get_old_sol(const unsigned int &local_dof_id,
                     const std::vector<unsigned int> &dof_indices_sys) {
    return d_sys.old_solution(dof_indices_sys[local_dof_id]);
  }

  /*!
   * @brief Reads old solution
   *
   * This function uses the supplied dof indices vector to compute solution
   *
   * @param local_dof_id Local dof id (eg 0 for constant element, {0,1,2,3,.
   * .} for linear element where numbers are local id of vertex)
   * @param local_var_id Local variable id (0 for tumor, 1 for chemical
   * potential)
   * @param dof_indices_sys_var Vector of Dof indices vector
   * @return solution Solution at the given dof
   */
  double get_old_sol_var(
    const unsigned int &local_dof_id, const unsigned int &local_var_id,
    const std::vector<std::vector<unsigned int>> &dof_indices_sys_var) {
    return d_sys.old_solution(dof_indices_sys_var[local_var_id][local_dof_id]);
  }

  /*!
   * @brief Set element solution in Libmesh system
   *
   * @param sol Solution vector
   * @param scale Factor if any
   */
  void set_elem_sol(const std::vector<double> &sol, double scale = 1.) {

    // Looping through elements
    for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

      init_dof(elem);
      d_sys.solution->set(get_global_dof_id(0), scale * sol[elem->id()]);
    }

    d_sys.solution->close();
    d_sys.update();
  }

  /*!
   * @brief Get element solution from Libmesh system
   *
   * @param sol Solution vector where we write the solution
   * @param scale Factor if any
   */
  void get_elem_sol(std::vector<double> &sol, double scale = 1.) {

    if (sol.size() < d_mesh.n_elem())
      sol = std::vector<double>(d_mesh.n_elem(), 0.);

    // Looping through elements
    for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

      init_dof(elem);
      sol[elem->id()] = scale * get_current_sol(0);
    }
  }

  /*!
   * @brief Initializes localized solution vector
   *
   * @param comm Reference to comm object
   */
  void init_localized_sol(const Parallel::Communicator &comm) {
    d_localized_sol = NumericVector<Number>::build(comm);
    d_localized_sol->init(d_sys.solution->size(), false, SERIAL);
    //if (comm.rank() == 0)
    //  d_localized_sol_std = std::vector<double>(d_sys.solution->size(), 0.);
  }

  /*!
   * @brief Localize the solution where the entries are number by the element id
   *
   * Here the fields are assumed to be constant in each element
   *
   * @param localize_sol Vector where solution are arranged element wise
   * @param var_ids Ids of variable for which we build the solutions
   * @param resize_vec Set to true to resize the vector
   */
  void localize_solution_with_elem_id_numbering_const_elem(
    std::vector<double> &localize_sol,
    std::vector<unsigned int> var_ids = {0}, bool resize_vec = true) {

    // TODO
    //  Should localize_to_one() be used so that only on processor 0 is localized
    //  Syntax: d_sys.solution->localize_to_one(d_localized_sol_std);

    // gather solution in all processors
    // sys.d_sys.current_local_solution->localize(collect_sol);
    d_sys.solution->localize(*d_localized_sol);

    if (d_mesh.comm().rank() > 0)
      return;

    // check if we need to resize the vector
    auto num_vars = var_ids.size();
    if (localize_sol.size() != d_mesh.n_elem() * num_vars) {
      if (resize_vec)
        localize_sol.resize(d_mesh.n_elem() * num_vars);
      else
        libmesh_error_msg(
          "localize_sol size should match collect_sol size for system " +
          d_sys_name);
    }

    // check if system has only 1 variable
    bool has_multiple_vars = d_num_vars > 1;
    if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
      libmesh_error_msg("Invalid var ids to collect and localize solution");

    for (const auto &elem : d_mesh.active_element_ptr_range()) {

      init_dof(elem);
      int counter = 0;
      for (auto var : var_ids) {

        // compute value at center
        auto val = 0.;
        if (has_multiple_vars)
          val += (*d_localized_sol)(get_var_global_dof_id(0, var));
        else
          val += (*d_localized_sol)(get_global_dof_id(0));

        localize_sol[elem->id() * num_vars + counter] = val;
        counter++;
      }
    }
  }

  /*!
   * @brief Localize the solution where the entries are number by the element id
   *
   * Here the fields are not piecewise constant
   *
   * @param localize_sol Vector where solution are arranged element wise
   * @param var_ids Ids of variable for which we build the solutions
   * @param resize_vec Set to true to resize the vector
   */
  void localize_solution_with_elem_id_numbering_non_const_elem(
    std::vector<double> &localize_sol,
    std::vector<unsigned int> var_ids = {0}, bool resize_vec = true) {

    // gather solution in all processors
    // sys.d_sys.current_local_solution->localize(collect_sol);
    d_sys.solution->localize(*d_localized_sol);

    if (d_mesh.comm().rank() > 0)
      return;

    // check if we need to resize the vector
    auto num_vars = var_ids.size();
    if (localize_sol.size() != d_mesh.n_elem() * num_vars) {
      if (resize_vec)
        localize_sol.resize(d_mesh.n_elem() * num_vars);
      else
        libmesh_error_msg(
          "localize_sol size should match collect_sol size for system " +
          d_sys_name);
    }

    // check if system has only 1 variable
    bool has_multiple_vars = d_num_vars > 1;
    if (!has_multiple_vars and (var_ids.size() > 1 or var_ids[0] != 0))
      libmesh_error_msg("Invalid var ids to collect and localize solution");

    for (const auto &elem : d_mesh.active_element_ptr_range()) {

      init_dof(elem);
      int counter = 0;
      for (auto var : var_ids) {

        // compute value at center
        auto val = 0.;
        for (unsigned int l = 0; l < d_phi.size(); l++) {
          if (has_multiple_vars)
            val += (*d_localized_sol)(get_var_global_dof_id(l, var));
          else
            val += (*d_localized_sol)(get_global_dof_id(l));
        }
        val = val / (double(d_phi.size()));

        localize_sol[elem->id() * num_vars + counter] = val;

        counter++;
      }
    }
  }

public:
  /*! @brief Name of system */
  std::string d_sys_name;

  /*! @brief Constant reference to system */
  TransientLinearImplicitSystem &d_sys;

  /*! @brief Mesh */
  const MeshBase &d_mesh;

  /*! @brief Number of variables in this system */
  unsigned int d_num_vars;

  /*! @brief Relative variable ids in the system for each variable */
  std::vector<unsigned int> d_var_id;

  /*! @brief Dof map of the system */
  const DofMap &d_dof_map_sys;

  /*! @brief Global dof ids at a given element */
  std::vector<unsigned int> d_dof_indices_sys;

  /*! @brief Global dof ids vector at a given element */
  std::vector<std::vector<unsigned int>> d_dof_indices_sys_var;

  /*! @brief Number of dofs in an element */
  unsigned int d_n_dofs;
  unsigned int d_n_dofs_var;

  /*! @brief FE type */
  FEType d_fe_type;

  /*! @brief Quadrature points information */
  QGauss d_qrule;

  /*! @brief FE object */
  UniquePtr<FEBase> d_fe;

  /*! @brief Reference to quadrature data */
  const std::vector<Point> &d_qpoints;
  const std::vector<Real> &d_JxW;
  const std::vector<std::vector<Real>> &d_phi;
  const std::vector<std::vector<RealGradient>> &d_dphi;

  /*! @brief Quadrature points information for integral over face */
  QGauss d_qrule_face;

  /*! @brief FE object */
  UniquePtr<FEBase> d_fe_face;

  /*! @brief Reference to quadrature data */
  const std::vector<Point> &d_qpoints_face;
  const std::vector<Real> &d_JxW_face;
  const std::vector<std::vector<Real>> &d_phi_face;
  const std::vector<std::vector<RealGradient>> &d_dphi_face;
  const std::vector<Point> &d_qface_normals;

  /*! @brief Local matrix and vector */
  DenseMatrix<Number> d_Ke;
  DenseVector<Number> d_Fe;

  /*! @brief Local matrix and vector for multi-variable systems */
  std::vector<std::vector<DenseSubMatrix<Number>>> d_Ke_var;
  std::vector<DenseSubVector<Number>> d_Fe_var;

  /*! @brief Localized solution vector */
  std::unique_ptr<NumericVector<Number>> d_localized_sol;
  std::vector<Number> d_localized_sol_std;

  /*! @brief init assembly matrix */
  bool d_assemble_matrix;

  /*! @brief Assembly type (implicit or explicit) */
  bool d_implicit_assembly;
};

} // namespace util

#endif // UTIL_SYSTEM_ABSTRACTION_H
