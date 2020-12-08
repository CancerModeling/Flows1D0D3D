////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NOISYCH_CAHNHILLIARD_H
#define NOISYCH_CAHNHILLIARD_H

#include "usystem/abstraction.hpp"
#include "usystem/stochastic_noise_assembly.hpp"

namespace noisych {

// forward declare
class Model;

/*! @brief Random initial condition for the concentration. */
Number initial_condition_cahnhilliard_random(const Point &p,
                                             const Parameters &es,
                                             const std::string &system_name,
                                             const std::string &var_name);

/*! @brief Smooth circular initial condition for the concentration. */
Number initial_condition_cahnhilliard_circle_2d(const Point &p,
                                                const Parameters &es,
                                                const std::string &system_name,
                                                const std::string &var_name);

Number initial_condition_cahnhilliard_circle_3d(const Point &p,
                                                const Parameters &es,
                                                const std::string &system_name,
                                                const std::string &var_name);

class CahnHilliardConfig {
public:
  unsigned int num_terms_bound;
  unsigned int seed;
  double scale;
  double length;
  double lower_bound;
  double upper_bound;

  double dt;
  double C_psi;
  double epsilon;
  double mobility_constant;

public:
  static CahnHilliardConfig from_parameter_file(const std::string &filename);

private:
  CahnHilliardConfig();
};

/*! @brief Class to perform assembly of the cahn-hilliard equation. */
class CahnHilliardAssembly : public util::BaseAssembly {

public:
  /*! @brief Constructor */
  CahnHilliardAssembly(
    const std::string &system_name,
    MeshBase &mesh,
    TransientLinearImplicitSystem &sys,
    const CahnHilliardConfig &config);

  /*! @brief Assembly function. Overrides the default assembly function */
  void assemble() override;

  /*! @brief Calculates new stochastic coefficients for our cylindrical Wiener process. */
  void calculate_new_stochastic_coefficients(double dt);

public:
  /*! @brief Assembles the noise from the cylindrical Wiener process. */
  util::StochasticNoiseAssembly d_noise_assembly;

  double d_dt;

  double d_C_psi;

  double d_epsilon;

  double d_mobility_constant;

private:
  /*! @brief Assembly */
  void assemble_1();
};

} // namespace noisych

#endif // NOISYCH_CAHNHILLIARD_H
