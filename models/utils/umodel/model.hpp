////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_BASEMODEL_H
#define UTIL_BASEMODEL_H

// Libmesh
#include "logger.hpp"
#include "qoi.hpp"
#include "utils.hpp"

// input deck
#include "uinp/inp.hpp"

// Assembly class
#include "usystem/abstraction.hpp"

// typedef input deck so that change its namespace or name does not effect
// the rest of the code
typedef util::InputDeck InpDeck;

namespace util {

/*!
 * @brief Abstraction of model class. Derived classes of this class can be
 * linked to network class easily with this abstraction.
 */
class BaseModel {

public:
  /*! @brief Constructor */
  BaseModel(Parallel::Communicator *comm, InpDeck &input, ReplicatedMesh &mesh,
            EquationSystems &tum_sys, util::Logger &log,
            const std::string &name)
      : d_name(name), d_step(0), d_time(0.), d_dt(0.), d_hmin(0.), d_hmax(0.),
        d_bounding_box(Point(), Point()), d_nonlinear_step(0),
        d_is_output_step(false), d_is_growth_step(false), d_log(log),
        d_delayed_msg(""), d_procRank(comm->rank()), d_procSize(comm->size()),
        d_input(input), d_comm_p(comm), d_mesh(mesh), d_tum_sys(tum_sys) {}

  /*! @brief Get equation system */
  const EquationSystems &get_system() const { return d_tum_sys; }
  EquationSystems &get_system() { return d_tum_sys; }

  const MeshBase &get_mesh() const { return d_mesh; }
  MeshBase &get_mesh() { return d_mesh; }

  /*! @brief Get input deck */
  const InpDeck &get_input_deck() const { return d_input; }
  InpDeck &get_input_deck() { return d_input; }

  /*! @brief Get mpi communicator */
  Parallel::Communicator *get_comm() const { return d_comm_p; }
  Parallel::Communicator *get_comm() { return d_comm_p; }

  /*! @brief Run model */
  virtual void run() = 0;

  /*! @brief Get various system classes */
  virtual BaseAssembly &get_assembly(const std::string &system) {
    libmesh_error_msg("Error: get_assembly should be defined in inheriting class");
  };

  virtual std::vector<util::BaseAssembly *> get_all_assembly() {
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
  std::pair<Point, Point> d_bounding_box;

  /*! @brief Current nonlinear step */
  unsigned int d_nonlinear_step;

  /*! @brief Is this output time step */
  bool d_is_output_step;

  /*! @brief Is this growth time step */
  bool d_is_growth_step;

  /*! @brief Reference to logger */
  util::Logger &d_log;

  /*! @brief List of systems */
  std::vector<std::string> d_sys_names;

  /*! @brief Store delayed message */
  std::string d_delayed_msg;

  /*! @brief Rank of processor */
  unsigned int d_procRank;

  /*! @brief Size of processors */
  unsigned int d_procSize;

protected:
  /*! @brief Output results of tumor system and network system */
  virtual void write_system(const unsigned int &t_step) = 0;

  /*! @brief Projects species concentrations to their physical range [0,1] */
  void project_solution_to_physical_range(const MeshBase &mesh,
                                          TransientLinearImplicitSystem &sys) {

    if (!d_input.d_project_solution_to_physical_range)
      return;

    // loop over dofs and project solution to [0,1] range
    const auto &sys_name = sys.name();

    unsigned int num_vars = 1;
    unsigned int sys_number = sys.number();
    unsigned int var_number = 0;
    if (sys_name == "Tumor")
      num_vars = 2;

    // loop over nodes and modify dofs
    for (const auto &node : mesh.local_node_ptr_range()) {

      const auto &dof = node->dof_number(sys_number, var_number, 0);

      auto val = sys.current_solution(dof);
      if (val < 0.)
        sys.solution->add(dof, -val);
      else if (val > 1.)
        sys.solution->add(dof, -val + 1.);
      else
        continue;
    }

    sys.solution->close();
    sys.update();
  }

  /*! @brief Solves tumor system */
  virtual void solve_system() = 0;

  /*! @brief Compute quantity of interest */
  virtual void compute_qoi() = 0;

  /*! @brief To store input parameters */
  InpDeck &d_input;

  /*! @brief Pointer to communicator */
  Parallel::Communicator *d_comm_p;

  /*! @brief Store the network mesh */
  ReplicatedMesh &d_mesh;

  /*! @brief Store the 2d/3d tumor system */
  EquationSystems &d_tum_sys;

  /*! @brief Compute and store QoI */
  QoIVec d_qoi;
};

} // namespace util

#endif // UTIL_BASEMODEL_H
