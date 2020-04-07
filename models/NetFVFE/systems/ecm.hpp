////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_ECM_H
#define NETFVFE_ECM_H

#include "abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*!
 * @brief Initial condition for ecm species
 *
 * @param p Point at which ic is computed
 * @param es Equation system
 * @param system_name Name of system
 * @param var_name Name of the variable
 * @param value Initial condition at given point
 */
Number initial_condition_ecm(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*!
 * @brief Class to perform assembly of ecm species
 */
class EcmAssembly : public BaseAssembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  EcmAssembly(Model * model, const std::string system_name, MeshBase &mesh,
      TransientLinearImplicitSystem & sys)
      : BaseAssembly(model, system_name, mesh, sys, 1,
                     {sys.variable_number("ecm")}) {}

  /*!
   * @brief Assembly function
   *
   * Overrides the default assembly function. It calls assemble_1,
   * assemble_2, or assemble_3 depending on user flag
   */
  void assemble() override;

private:

  /*!
   * @brief Assembly over volume of element
   *
   * If assembly method == 1
   * In this we simply implement the assembly under iterative nonlinear
   * scheme. In source terms, those which are linear with respect to system
   * variable, we consider implicit scheme.
   *
   * If assembly method == 2
   * Same as above, but now we project the species concentration to
   * physical range [0,1].
   */
  void assemble_1();
};

} // namespace netfvfe

#endif // NETFVFE_ECM_H
