////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_GRAD_TAF_H
#define NETFV_GRAD_TAF_H

#include "abstraction.hpp"

namespace netfv {

// forward declare
class Model;

/*!
 * @brief Class to perform assembly of gradient of TAF
 */
class GradTafAssembly : public BaseAssembly {

public:
  /*!
   * @brief Constructor
   *
   * @param model Model class
   * @param sys_name Name of system
   * @param sys System
   */
  GradTafAssembly(Model * model, const std::string system_name,
      TransientLinearImplicitSystem & sys)
      : BaseAssembly(model, system_name, sys, 2,
                     {sys.variable_number("taf_gradx"),
                      sys.variable_number("taf_grady")}) {

    if (sys.get_mesh().mesh_dimension() > 2) {
      d_num_vars = 3;
      d_var_id.push_back(d_sys.variable_number("taf_gradx"));
      d_dof_indices_sys_var.resize(d_num_vars);
    }
  }

  /*!
   * @brief Assembly function
   *
   * Overrides the default assembly function.
   *
   * For grad taf, we do not need to solve. We have another function which
   * computes the grad of taf at center of elements and updates the solution
   * data.
   */
  void assemble() override {}

  /*!
   * @brief Solve for gradient taf
   */
  void solve();
};

} // namespace netfv

#endif // NETFV_GRAD_TAF_H
