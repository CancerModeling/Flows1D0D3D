////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef AVAFV_ABSTRACTION_H
#define AVAFV_ABSTRACTION_H

#include "utilLibs.hpp"
#include "utils.hpp"

#include <string>

namespace avafv {

// forward declare
class Model;

/*!
 * @brief Class to perform assembly of tumor species
 */
class BaseAssembly : public System::Assembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  BaseAssembly(Model *model, const std::string system_name,
               TransientLinearImplicitSystem &sys, const unsigned int &num_vars,
               const std::vector<unsigned int> &var_id)
      : d_sys_name(system_name), d_model_p(model), d_sys(sys),
        d_num_vars(num_vars), d_var_id(var_id),
        d_dof_map_sys(d_sys.get_dof_map()) {

    d_dof_indices_sys_var.resize(d_num_vars);
  }

  /*!
   * @brief Assembly function
   *
   * To be overridden by inheriting classes
   */
  virtual void assemble() = 0;

  /*!
   * @brief Initializes the dof state with given element
   *
   * @param elem Pointer to the element
   */
  void init_dof(const Elem * elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * For system with more than one variable
   *
   * @param elem Pointer to the element
   */
  void init_var_dof(const Elem * elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
    for (unsigned int var = 0; var < d_num_vars; var++)
      d_dof_map_sys.dof_indices(elem, d_dof_indices_sys_var[var], d_var_id[var]);
  }

  /*!
   * @brief Initializes the dof state with given element
   *
   * It populates the dof_indices passed to the function
   *
   * @param elem Pointer to the element
   * @param dof_indices_sys Dof indices vector
   */
  void init_dof(const Elem * elem, std::vector<unsigned int> &dof_indices_sys) {

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
  void init_var_dof(const Elem *elem, std::vector<unsigned int>
      &dof_indices_sys,
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
    return d_sys.current_solution(d_dof_indices_sys_var[local_var_id][local_dof_id]);
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
  double get_current_sol_var(const unsigned int &local_dof_id,
                         const unsigned int &local_var_id,
                         const std::vector<std::vector<unsigned int>> &dof_indices_sys_var) {
    return d_sys.current_solution(dof_indices_sys_var[local_var_id][local_dof_id]);
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
    return d_sys.old_solution(d_dof_indices_sys_var[local_var_id][local_dof_id]);
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
  double get_old_sol_var(const unsigned int &local_dof_id,
                         const unsigned int &local_var_id,
                         const std::vector<std::vector<unsigned int>> &dof_indices_sys_var) {
    return d_sys.old_solution(dof_indices_sys_var[local_var_id][local_dof_id]);
  }


public:

  /*! @brief Name of system */
  std::string d_sys_name;

  /*! @brief Pointer reference to model */
  Model *d_model_p;

  /*! @brief Constant reference to system */
  TransientLinearImplicitSystem &d_sys;

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
};

} // namespace avafv

#endif // AVAFV_ABSTRACTION_H
