////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_ABSTRACTION_H
#define NETFVFE_ABSTRACTION_H

#include "utilLibs.hpp"
#include "utils.hpp"
#include "../inp/inp.hpp"

#include <string>

namespace netfvfe {

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
  BaseAssembly(Model *model, const std::string system_name, MeshBase &mesh,
               TransientLinearImplicitSystem &sys, const unsigned int &num_vars,
               const std::vector<unsigned int> &var_id)
      : d_sys_name(system_name), d_model_p(model), d_mesh(mesh), d_sys(sys),
        d_num_vars(num_vars), d_var_id(var_id),
        d_n_dofs(0), d_n_dofs_var(0),
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
        d_qface_normals(d_fe_face->get_normals()) {

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
   * @brief Initializes the fe and local matrix and vector
   *
   * @param elem Pointer to the element
   */
  void init_fe(const Elem * elem) {

    d_fe->reinit(elem);

    d_Fe.resize(d_n_dofs);
    d_Ke.resize(d_n_dofs, d_n_dofs);


    if (d_num_vars > 1)
      for (unsigned int i=0; i<d_num_vars; i++) {

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
    void init_face_fe(const Elem * elem, unsigned short side) {
        d_fe_face->reinit(elem, side);
    }

  /*!
   * @brief Initializes the dof state with given element
   *
   * @param elem Pointer to the element
   */
  void init_dof(const Elem * elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
    d_n_dofs = d_dof_indices_sys.size();
    d_n_dofs_var = d_n_dofs;

    if (d_num_vars > 1) {
      for (unsigned int var = 0; var < d_num_vars; var++)
        d_dof_map_sys.dof_indices(elem, d_dof_indices_sys_var[var], d_var_id[var]);

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
  void init_var_dof(const Elem * elem) {

    d_dof_map_sys.dof_indices(elem, d_dof_indices_sys);
    for (unsigned int var = 0; var < d_num_vars; var++)
      d_dof_map_sys.dof_indices(elem, d_dof_indices_sys_var[var], d_var_id[var]);

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
};

} // namespace netfvfe

#endif // NETFVFE_ABSTRACTION_H
