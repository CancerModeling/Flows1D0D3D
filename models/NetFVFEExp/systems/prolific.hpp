////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFEEXP_PROLIFIC_H
#define NETFVFEEXP_PROLIFIC_H

#include "usystem/abstraction.hpp"

namespace netfvfeexp {

// forward declare
class Model;

/*! @brief Initial condition for tumor species */
Number initial_condition_pro(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of tumor species */
class ProAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  ProAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 2,
                           {sys.variable_number("prolific"),
                            sys.variable_number("chemical_prolific")}),
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

} // namespace netfvfeexp

#endif // NETFVFEEXP_PROLIFIC_H
