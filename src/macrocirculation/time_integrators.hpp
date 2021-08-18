////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TIME_INTEGRATORS_H
#define TUMORMODELS_TIME_INTEGRATORS_H

#include <vector>

namespace macrocirculation {

// forward declarations:
class RightHandSideEvaluator;

struct ButcherScheme {
  /*! @brief The lower diagonal of the A-matrix of the Butcher scheme. */
  std::vector<double> a;
  /*! @brief The b vector of the Butcher scheme. */
  std::vector<double> b;
  /*! @brief The c vector of the Butcher scheme. */
  std::vector<double> c;
};

/*! @brief Creates the butcher scheme for the explicit euler. */
ButcherScheme create_explicit_euler();

/*! @brief Creates the butcher scheme for the 3rd order ssp method. */
ButcherScheme create_ssp_method();

/*! @brief Checks the consistency of a butcher scheme. */
bool check_consistency_rkm(const ButcherScheme &bs);

/*! @brief Class for evaluating butcher schemes on gmm vectors. */
class TimeIntegrator {
public:
  TimeIntegrator(ButcherScheme bs, std::size_t num_dofs);

  void apply(const std::vector<double> &u_prev, double t, double tau, RightHandSideEvaluator &rhs, std::vector<double> &u_now) const;

private:
  ButcherScheme d_bs;

  mutable std::vector<std::vector<double>> d_k;
  mutable std::vector<double> d_tmp;
};

} // namespace macrocirculation

#endif //TUMORMODELS_TIME_INTEGRATORS_H
