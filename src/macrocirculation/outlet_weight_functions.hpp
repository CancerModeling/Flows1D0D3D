////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_OUTLET_WEIGHT_FUNCTIONS_HPP
#define TUMORMODELS_OUTLET_WEIGHT_FUNCTIONS_HPP

#include "libmesh/point.h"
#include <vector>

namespace lm = libMesh;

namespace macrocirculation {

/*! @brief Base class to define outlet weight functions */
class BaseOutletRadial {
public:
  BaseOutletRadial(int type, lm::Point x, double r) : d_type(type), d_x(x), d_r(r), d_c(1.) {}

  /* @brief set the normalizing constant */
  void set_normalize_const(double c) {d_c = c; }

  /* @brief type of radial function. 0 - const, 1 - linear, 2 - quadratic */
  int d_type;

  /* @brief outlet point */
  lm::Point d_x;

  /* @brief radius of ball on which this function is supported */
  double d_r;

  /* @brief normalizing constant */
  double d_c;

  virtual double operator()(const lm::Point &x) const {
    std::cerr << "Error: Function should be defined by inheriting class.\n";
    exit(EXIT_FAILURE);
  }
};

/*! @brief Constant outlet weight functions - f(r) = 1, r in [0,1] */
class ConstOutletRadial : public BaseOutletRadial {
public:
  ConstOutletRadial(lm::Point x, double r) : BaseOutletRadial(0, x, r) {
    //d_c = 3. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c;
    else
      return 0.;
  }
};

/*! @brief Linear outlet weight functions - f(r) = 1 - r, r in [0,1] */
class LinearOutletRadial : public BaseOutletRadial {
public:
  LinearOutletRadial(lm::Point x, double r) : BaseOutletRadial(0, x, r) {
    //d_c = 12. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c * (1. - r);
    else
      return 0.;
  }
};

/*! @brief Gaussian outlet weight functions - f(r) = exp[-r^2/(2*sigma^2)], r in [0,1] */
class GaussianOutletRadial : public BaseOutletRadial {
public:
  GaussianOutletRadial(lm::Point x, double r, double sigma) : BaseOutletRadial(0, x, r), d_sigma(sigma) {
    //double m2 = 0.5 * d_sigma * d_sigma * (d_sigma * std::sqrt(M_PI * 2.) * std::erf(1. / (d_sigma * std::sqrt(2))) - 2. * std::exp(-1. / (2. * d_sigma * d_sigma)));
    //d_c = 1. / (4. * M_PI * std::pow(d_r, 3) * m2); // normalizing constant
  }

  /* @brief std of gaussian */
  double d_sigma;

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c * std::exp(-std::pow(r, 2) / (2. * d_sigma * d_sigma));
    else
      return 0.;
  }
};


} // namespace macrocirculation

#endif //TUMORMODELS_OUTLET_WEIGHT_FUNCTIONS_HPP
