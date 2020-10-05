////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef AVAFV_GRAD_TAF_H
#define AVAFV_GRAD_TAF_H

#include "usystem/abstraction.hpp"

namespace avafv {

// forward declare
class Model;

/*! @brief Class to perform assembly of gradient of TAF */
class GradTafAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  GradTafAssembly(Model *model, const std::string system_name, MeshBase &mesh,
                  TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(
          system_name, mesh, sys, mesh.mesh_dimension(),
          (mesh.mesh_dimension() == 2
             ? std::vector<unsigned int>{sys.variable_number("taf_gradx"),
                                         sys.variable_number("taf_grady")}
             : std::vector<unsigned int>{sys.variable_number("taf_gradx"),
                                         sys.variable_number("taf_grady"),
                                         sys.variable_number(
                                           "taf_gradz")})),
        d_model_p(model) {}

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override {}

  /*! @brief Directly computes the gradient of taf and sets the values of
   * solution */
  void solve();

public:
  /*! @brief Pointer reference to model */
  Model *d_model_p;
};

} // namespace avafv

#endif // AVAFV_GRAD_TAF_H
