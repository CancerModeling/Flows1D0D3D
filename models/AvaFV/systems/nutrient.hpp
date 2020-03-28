////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef AVAFV_NUTRIENT_H
#define AVAFV_NUTRIENT_H

#include "abstraction.hpp"

namespace avafv {

// forward declare
class Model;

/*!
 * @brief Initial condition for nutrient
 *
 * @param p Point at which ic is computed
 * @param es Equation system
 * @param system_name Name of system
 * @param var_name Name of the variable
 * @param value Initial condition at given point
 */
Number initial_condition_nut(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*!
 * @brief Boundary condition for nutrient
 *
 * @param es Equation system
 */
void boundary_condition_nut(EquationSystems &es);

/*!
 * @brief Class to perform assembly of pressure in tissue domain
 */
class NutAssembly : public BaseAssembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  NutAssembly(Model * model, const std::string system_name,
      TransientLinearImplicitSystem & sys)
      : BaseAssembly(model, system_name, sys, 1,
                     {sys.variable_number("nutrient")}) {}

  /*!
   * @brief Assembly function
   *
   * Overrides the default assembly function. It calls assemble_1,
   * assemble_2, or assemble_3 depending on user flag
   */
  void assemble() override;

private:

  /*!
   * @brief Assembly function 2
   *
   * Same as assemble_1, but now we project the species concentration to
   * physical range [0,1].
   */
  void assemble_vol();

  /*!
   * @brief Assembly of terms over face of element
   */
  void assemble_face();
};

} // namespace avafv

#endif // AVAFV_NUTRIENT_H
