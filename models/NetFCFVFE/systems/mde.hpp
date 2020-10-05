////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFCFVFE_MDE_H
#define NETFCFVFE_MDE_H

#include "usystem/abstraction.hpp"

namespace netfcfvfe {

// forward declare
class Model;

/*! @brief Initial condition for MDE */
Number initial_condition_mde(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of MDE */
class MdeAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  MdeAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 1,
                           {sys.variable_number("mde")}),
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

} // namespace netfcfvfe

#endif // NETFCFVFE_MDE_H
