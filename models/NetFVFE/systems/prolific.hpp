////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_PROLIFIC_H
#define NETFVFE_PROLIFIC_H

#include "usystem/abstraction.hpp"
#include "usystem/stochastic_noise_assembly.hpp"

namespace netfvfe {

// forward declare
class Model;

/*! @brief Initial condition for tumor species */
Number initial_condition_pro(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/*! @brief Class to perform assembly of tumor species */
class ProAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  ProAssembly(Model *model,
              const std::string &system_name,
              MeshBase &mesh,
              TransientLinearImplicitSystem &sys);

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override;

  /*! @brief Calculates new stochastic coefficients for our cylindrical Wiener process. */
  void calculate_new_stochastic_coefficients(double dt);

public:
  /*! @brief Pointer reference to model */
  Model *d_model_p;

  /*! @brief Assembles the noise from the cylindrical Wiener process. */
  util::StochasticNoiseAssembly d_noise_assembly;

private:
  /*! @brief Assembly */
  void assemble_1();
};

} // namespace netfvfe

#endif // NETFVFE_PROLIFIC_H
