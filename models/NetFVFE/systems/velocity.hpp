////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_VELOCITY_H
#define NETFVFE_VELOCITY_H

#include "usystem/abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Class to perform assembly of velocity */
class VelAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  VelAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(
          system_name, mesh, sys, mesh.mesh_dimension(),
          (mesh.mesh_dimension() == 2
             ? std::vector<unsigned int>{sys.variable_number("velocity_x"),
                                         sys.variable_number("velocity_y")}
             : std::vector<unsigned int>{sys.variable_number("velocity_x"),
                                         sys.variable_number("velocity_y"),
                                         sys.variable_number(
                                           "velocity_z")})),
        d_model_p(model) {}

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override;

public:
  /*! @brief Pointer reference to model */
  Model *d_model_p;

private:
  /*! @brief Assembly */
  void assemble_1();
};

} // namespace netfvfe

#endif // NETFVFE_VELOCITY_H
