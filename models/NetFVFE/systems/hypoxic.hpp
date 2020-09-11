////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_HYPOXIC_H
#define NETFVFE_HYPOXIC_H

#include "usystem/abstraction.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Initial condition for hypoxic species */
Number initial_condition_hyp(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of hypoxic species */
class HypAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  HypAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 2,
                           {sys.variable_number("hypoxic"),
                            sys.variable_number("chemical_hypoxic")}),
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
