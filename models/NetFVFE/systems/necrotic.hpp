////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_NECROTIC_H
#define NETFVFE_NECROTIC_H

#include "usystem/abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Initial condition for necrotic species */
Number initial_condition_nec(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of necrotic species */
class NecAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  NecAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 1,
                           {sys.variable_number("necrotic")}),
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

#endif // NETFVFE_HYPOXIC_H
