////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_OUTPUT_INTERPOLATOR_HPP
#define TUMORMODELS_OUTPUT_INTERPOLATOR_HPP

#include <vector>

namespace macrocirculation {

class Interpolator {
public:
  /** @brief Adds a value at a certain time.  */
  void add(double ptime, double pvalue);

  /** @brief Returns the inter- or extrapolated value at the given time. We use the previous values provided with add for the interpolation.  */
  double operator()(double current_time) const;

private:
  std::vector<double> time;
  std::vector<double> value;
};

/** @brief Interpolates a quantity of interest at several consecutively numbered points called outlets.  */
class OutletInterpolator {
public:
  /** @brief Creates an Interpolater for several outlets. */
  OutletInterpolator(int num_outlets);

  /** @brief Adds a value at one outlet at a certain time.  */
  void add(int outlet, double time, double value);

  /** @brief Returns the inter- or extrapolated value at the outlet at the given time. We use the the previous values provided with add for the interpolation. */
  double operator()(int outlet, double current_time) const;

private:
  std::vector<Interpolator> value_interpolator;
};

} // namespace macrocirculation

#endif