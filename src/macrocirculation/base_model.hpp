////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "libmesh_includes.hpp"
#include "logger.hpp"
#include "assembly_system.hpp"

namespace macrocirculation {

/*!
 * @brief Abstraction of model class. Derived classes of this class can be
 * linked to network class easily with this abstraction.
 */
class BaseModel {

public:
  /*! @brief Constructor */
  BaseModel(lm::Parallel::Communicator *comm, lm::ReplicatedMesh &mesh,
            lm::EquationSystems &eq_sys, Logger &log,
            const std::string &name)
      : d_name(name), d_step(0), d_time(0.), d_dt(0.), d_hmin(0.), d_hmax(0.),
        d_bounding_box(lm::Point(), lm::Point()), 
        d_nonlinear_step(0), d_is_output_step(false), 
        d_log(log), d_procRank(comm->rank()), d_procSize(comm->size()),
        d_comm_p(comm), d_mesh(mesh), d_eq_sys(eq_sys) {}

  /*! @brief Get equation system */
  const lm::EquationSystems &get_system() const { return d_eq_sys; }
  lm::EquationSystems &get_system() { return d_eq_sys; }

  const lm::MeshBase &get_mesh() const { return d_mesh; }
  lm::MeshBase &get_mesh() { return d_mesh; }

  /*! @brief Get mpi communicator */
  lm::Parallel::Communicator *get_comm() const { return d_comm_p; }
  lm::Parallel::Communicator *get_comm() { return d_comm_p; }

  /*! @brief Get various system classes */
  virtual BaseAssembly &get_assembly(const std::string &system) {
    libmesh_error_msg("Error: get_assembly should be defined in inheriting class");
  };

  virtual std::vector<BaseAssembly *> get_all_assembly() {
    libmesh_error_msg("Error: get_all_assembly should be defined in inheriting "
                      "class");
  }

public:
  /*! @brief To store input parameters */
  std::string d_name;

  /*! @brief To store input parameters */
  unsigned int d_step;

  /*! @brief Current time */
  double d_time;

  /*! @brief Current time step */
  double d_dt;

  /*! @brief hmax and hmin */
  double d_hmin;
  double d_hmax;

  /*! @brief Bounding box */
  std::pair<lm::Point, lm::Point> d_bounding_box;

  /*! @brief Current nonlinear step */
  unsigned int d_nonlinear_step;

  /*! @brief Is this output time step */
  bool d_is_output_step;

  /*! @brief Is this growth time step */
  bool d_is_growth_step;

  /*! @brief Reference to logger */
  Logger &d_log;

  /*! @brief List of systems */
  std::vector<std::string> d_sys_names;

  /*! @brief Rank of processor */
  unsigned int d_procRank;

  /*! @brief Size of processors */
  unsigned int d_procSize;

protected:

  /*! @brief Pointer to communicator */
  lm::Parallel::Communicator *d_comm_p;

  /*! @brief Store the network mesh */
  lm::ReplicatedMesh &d_mesh;

  /*! @brief Store the 2d/3d tumor system */
  lm::EquationSystems &d_eq_sys;
};

} // namespace macrocirculation

#endif // BASE_MODEL_H
