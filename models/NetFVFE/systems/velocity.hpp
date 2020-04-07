////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_VELOCITY_H
#define NETFVFE_VELOCITY_H

#include "abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*!
 * @brief Class to perform assembly of velocity
 */
class VelAssembly : public BaseAssembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  VelAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : BaseAssembly(
            model, system_name, mesh, sys, mesh.mesh_dimension(),
            (mesh.mesh_dimension() == 2
                 ? std::vector<unsigned int>{sys.variable_number("velocity_x"),
                                             sys.variable_number("velocity_y")}
                 : std::vector<unsigned int>{
                       sys.variable_number("velocity_x"),
                       sys.variable_number("velocity_y"),
                       sys.variable_number("velocity_z")})) {}

  /*!
   * @brief Assembly function
   *
   * Overrides the default assembly function.
   */
  void assemble() override;
};

} // namespace netfvfe

#endif // NETFVFE_VELOCITY_H
