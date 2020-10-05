////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_TAF_H
#define NETFVFE_TAF_H

#include "usystem/abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Initial condition for TAF */
Number initial_condition_taf(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of TAF */
class TafAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  TafAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 1,
                           {sys.variable_number("taf")}),
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

#endif // NETFVFE_TAF_H
