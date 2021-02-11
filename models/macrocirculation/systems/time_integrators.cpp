////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "time_integrators.hpp"
#include "right_hand_side_evaluator.hpp"

#include "gmm.h"

namespace macrocirculation {

ButcherScheme create_explicit_euler() {
  ButcherScheme bs;
  bs.c.push_back(0);
  bs.b.push_back(1);
  return bs;
}

ButcherScheme create_ssp_method() {
  ButcherScheme bs;
  // c vector:
  bs.c.push_back(0);
  bs.c.push_back(1);
  bs.c.push_back(0.5);
  // b vector:
  bs.b.push_back(1. / 6.);
  bs.b.push_back(1. / 6.);
  bs.b.push_back(4. / 6.);
  // a matrix:
  bs.a.push_back(1.);
  bs.a.push_back(0.25);
  bs.a.push_back(0.25);
  return bs;
}

bool check_consistency_rkm(const ButcherScheme &bs) {
  throw std::runtime_error("not implemented yet.");
}

TimeIntegrator::TimeIntegrator(ButcherScheme bs, std::size_t num_dofs)
    : d_bs(std::move(bs)),
      d_k(std::vector<std::vector<double>>(bs.b.size(), std::vector<double>(num_dofs, 0))),
      d_tmp(num_dofs, 0) {}

template<std::size_t degree>
void TimeIntegrator::apply(const std::vector<double> &u_prev,
                           const double t,
                           const double tau,
                           RightHandSideEvaluator<degree> &rhs,
                           std::vector<double> &u_now) const {
  // evaluate ks
  for (std::size_t i = 0; i < d_bs.b.size(); i += 1) {
    double t_s = t + tau * d_bs.c[i];
    d_tmp = u_prev;
    for (std::size_t j = 0; j < i; j += 1) {
      std::size_t coeff = tau * d_bs.a.at(i * (i - 1) + j);
      gmm::add(d_tmp, gmm::scaled(d_k[j], coeff), d_tmp);
    }
    rhs.evaluate(t_s, d_tmp, d_k[i]);
  }
  // evaluate bs
  u_now = u_prev;
  for (std::size_t i = 0; i < d_bs.b.size(); i += 1) {
    std::size_t coeff = tau * d_bs.b.at(i);
    gmm::add(u_now, gmm::scaled(d_k[i], coeff), u_now);
  }
}

template void TimeIntegrator::apply<0>(const std::vector<double> &, const double, const double, RightHandSideEvaluator<0> &, std::vector<double> &) const;
template void TimeIntegrator::apply<1>(const std::vector<double> &, const double, const double, RightHandSideEvaluator<1> &, std::vector<double> &) const;
template void TimeIntegrator::apply<2>(const std::vector<double> &, const double, const double, RightHandSideEvaluator<2> &, std::vector<double> &) const;
template void TimeIntegrator::apply<3>(const std::vector<double> &, const double, const double, RightHandSideEvaluator<3> &, std::vector<double> &) const;

} // namespace macrocirculation