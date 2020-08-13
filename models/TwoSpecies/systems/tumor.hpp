////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TWOSP_TUMOR_H
#define TWOSP_TUMOR_H

#include "usystem/abstraction.hpp"

namespace twosp {

// forward declare
class Model;

/*! @brief Initial condition for tumor species */
Number initial_condition_tum(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of tumor species */
class TumAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  TumAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 2,
                           {sys.variable_number("tumor"),
                            sys.variable_number("chemical_tumor")}),
        d_model_p(model) {}

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override ;

public:

  /*! @brief Pointer reference to model */
  Model *d_model_p;

private:

  /*! @brief Assembly */
  void assemble_1();
};

} // namespace twosp

#endif // TWOSP_TUMOR_H
