////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2022 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "outlet_interpolator.hpp"
#include <cassert>
#include <stdexcept>

namespace macrocirculation {

void Interpolator::add(double ptime, double pvalue) {
  time.push_back(ptime);
  value.push_back(pvalue);
}

double Interpolator::operator()(double current_time) const {
  if (value.size() < 1)
    throw std::runtime_error("Not enough values to interpolate");

  double value_now = value[value.size() - 1];
  double value_prev = value[value.size() - 2];
  double time_now = time[time.size() - 1];
  double time_prev = time[time.size() - 2];
  double zeta = (current_time - time_prev) / (time_now - time_prev);
  double value_interp = value_prev + zeta * (value_now - value_prev);
  return value_interp;
}


OutletInterpolator::OutletInterpolator(int num_outlets)
    : value_interpolator(num_outlets) {}

void OutletInterpolator::add(int outlet, double time, double value) {
  assert(outlet < value_interpolator.size());
  value_interpolator[outlet].add(time, value);
}

double OutletInterpolator::operator()(int outlet, double current_time) const {
  assert(outlet < value_interpolator.size());
  return value_interpolator[outlet](current_time);
}

} // namespace macrocirculation