////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_PRESSURE_H
#define NETFVFE_PRESSURE_H

#include "usystem/abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Initial condition for pressure */
Number initial_condition_pres(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Boundary condition for pressure */
void boundary_condition_pres(EquationSystems &es);

/*! @brief Class to perform assembly of pressure in tissue domain */
class PressureAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  PressureAssembly(Model *model, const std::string system_name, MeshBase &mesh,
                   TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 1,
                     {sys.variable_number("pressure")}), d_model_p(model) {}

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override;

public:

  /*! @brief Pointer reference to model */
  Model *d_model_p;

private:

  /*! @brief Assembly */
  void assemble_1();

  /*! @brief Assemble coupling between 3d and 1d pressure */
  void assemble_1d_coupling();

  /*! @brief Assembly of terms over face of element */
  void assemble_face();
};

} // namespace netfvfe

#endif // NETFVFE_PRESSURE_H
