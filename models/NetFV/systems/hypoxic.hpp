////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_HYPOXIC_H
#define NETFV_HYPOXIC_H

#include "usystem/abstraction.hpp"

namespace netfv {

// forward declare
class Model;

/*! @brief Initial condition for hypoxic species */
Number initial_condition_hyp(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief A kernel function for initial condition */
Number initial_condition_hyp_kernel(const Point &p,
                                    const unsigned int &dim,
                                    const std::string &ic_type,
                                    const std::vector<double> &ic_center,
                                    const std::vector<double> &tum_ic_radius,
                                    const std::vector<double> &hyp_ic_radius);

/*! @brief Class to perform assembly of hypoxic species */
class HypAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  HypAssembly(Model *model, const std::string system_name, MeshBase &mesh,
              TransientLinearImplicitSystem &sys)
      : util::BaseAssembly(system_name, mesh, sys, 1,
                           {sys.variable_number("hypoxic")}),
        d_model_p(model) {}

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override;

public:
  /*! @brief Pointer reference to model */
  Model *d_model_p;

private:
  /*! @brief Assembly */
  void assemble_1();

  /*! @brief Assembly of terms over face of element */
  void assemble_face();
};

} // namespace netfv

#endif // NETFV_HYPOXIC_H
