////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NOISYCH_MODEL_H
#define NOISYCH_MODEL_H

#include "modelUtil.hpp"
#include "systems/systems.hpp"

namespace noisych {

void model_setup_run( const std::string &filename, Parallel::Communicator *comm);

/*!
 * @brief Simple model for a simplified Cahn-Hilliard equation.
 */
class Model {

public:
  /*! @brief Constructor */
  Model(const std::string &filename,
        Parallel::Communicator *comm,
        Mesh &mesh,
        EquationSystems &sys,
        TransientLinearImplicitSystem &ch);

  /*! @brief Run model. */
  void run();

private:
  void write_system(const unsigned int &t_step);

private:
  /*! @brief The equation system. */
  EquationSystems& d_sys;

  /*! @brief Current mesh. */
  Mesh& d_mesh;

  /*! @brief Assembly objects. */
  CahnHilliardAssembly d_ch;

  /*! @brief Current time step. */
  std::size_t d_step;

  /*! @brief Total number of time steps. */
  std::size_t d_steps;

  /*! @brief Current time. */
  double d_time;

  /*! @brief Current time step size. */
  double d_dt;
};

} // namespace noisych

#endif // NOISYCH_MODEL_H
