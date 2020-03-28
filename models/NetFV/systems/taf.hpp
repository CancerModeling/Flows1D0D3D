////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_TAF_H
#define NETFV_TAF_H

#include "abstraction.hpp"

namespace netfv {

// forward declare
class Model;

/*!
 * @brief Initial condition for TAF
 *
 * @param p Point at which ic is computed
 * @param es Equation system
 * @param system_name Name of system
 * @param var_name Name of the variable
 * @param value Initial condition at given point
 */
Number initial_condition_taf(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*!
 * @brief Class to perform assembly of TAF
 */
class TafAssembly : public BaseAssembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  TafAssembly(Model * model, const std::string system_name,
      TransientLinearImplicitSystem & sys)
      : BaseAssembly(model, system_name, sys, 1,
                     {sys.variable_number("taf")}) {}

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
   * In this we simply implement the assembly under iterative nonlinear
   * scheme. In source terms, those which are linear with respect to system
   * variable, we consider implicit scheme.
   */
  void assemble_1();

  /*!
   * @brief Assembly over volume of element
   *
   * Same as assemble_1, but now we project the species concentration to
   * physical range [0,1].
   */
  void assemble_2();

  /*!
   * @brief Assembly over volume of element
   *
   * Same as assemble_2, but now all the terms in source term are handled
   * explicitly.
   */
  void assemble_3();

  /*!
   * @brief Assembly of terms over face of element
   */
  void assemble_face();
};

} // namespace netfv

#endif // NETFV_TAF_H
