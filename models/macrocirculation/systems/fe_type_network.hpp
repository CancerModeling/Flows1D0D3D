////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_FE_TYPE_NETWORK_H
#define TUMORMODELS_FE_TYPE_NETWORK_H

#include <utility>

#include "libmesh/libmesh.h"
#include "graph_storage.hpp"

namespace macrocirculation {

constexpr double legendre2(double x) { return 0.5 * (3 * x * x - 1); }

struct QuadratureFormula {
  std::vector<double> ref_points;
  std::vector<double> ref_weights;

  std::size_t size() const { return ref_points.size(); }
};

inline QuadratureFormula create_gauss3() {
  QuadratureFormula qf;
  qf.ref_points = {-0.774596669241, 0, 0.774596669241};
  qf.ref_weights = {+5. / 9., 8. / 9., +5. / 9.};
  return qf;
}

inline QuadratureFormula create_trapezoidal_rule() {
  QuadratureFormula qf;
  qf.ref_points = {-1, +1};
  qf.ref_weights = {1, 1};
  return qf;
}

inline QuadratureFormula create_midpoint_rule() {
  QuadratureFormula qf;
  qf.ref_points = {0};
  qf.ref_weights = {2.};
  return qf;
}

/*
class FETypeNetwork {
public:
  explicit FETypeNetwork(QuadratureFormula qf)
      : d_qf(std::move(qf)),
        d_phi(3, std::vector<double>(d_qf.size())),
        d_dphi(3, std::vector<double>(d_qf.size())),
        d_JxW(d_qf.size()) {
    for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
      d_phi[0][qp] = 1;
      d_phi[1][qp] = d_qf.ref_points[qp];
      d_phi[2][qp] = legendre2(d_qf.ref_points[qp]);

      d_dphi[0][qp] = 0;
      d_dphi[1][qp] = NAN;
      d_dphi[2][qp] = NAN;

      d_JxW[qp] = NAN;
    }
  }

  void reinit(const Edge &e) {
    auto length = e.get_length();

    for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
      d_dphi[1][qp] = length / 2;
      d_dphi[2][qp] = 6 * d_qf.ref_points[qp] / length;
      d_JxW[qp] = d_qf.ref_weights[qp] * length / 2.;
    }
  };

  const std::vector<std::vector<double>> &get_phi() const { return d_phi; };

  const std::vector<std::vector<double>> &get_dphi() const { return d_dphi; };

  const std::vector<double> &get_JxW() const { return d_JxW; };

private:
  QuadratureFormula d_qf;

  std::vector<std::vector<double>> d_phi;

  std::vector<std::vector<double>> d_dphi;

  std::vector<double> d_JxW;
};

class FETypeInnerBdryNetwork {
public:
  explicit FETypeInnerBdryNetwork()
    : d_phi_r(3), d_phi_l(3) {
    d_phi_l[0] = 1;
    d_phi_r[0] = 1;

    d_phi_l[1] = NAN;
    d_phi_r[1] = NAN;

    d_phi_l[2] = NAN;
    d_phi_r[2] = NAN;
  }

  /// Assume: ----e_left----(v)----e_right----
  void reinit(const Vertex & v, const Edge &e_left, const Edge &e_right) {
    const double orientation_left = e_left.get_vertex_neighbors()[1] == v.get_id() ? +1 : -1;
    const double orientation_right = e_right.get_vertex_neighbors()[0] == v.get_id() ? +1 : -1;

    const double point_left = +1*orientation_left;
    const double point_right = -1*orientation_right;

    d_phi_l[1] = point_left;
    d_phi_r[1] = point_right;

    d_phi_l[2] = legendre2(point_left);
    d_phi_r[2] = legendre2(point_right);
  };

  const std::vector<double> &get_phi_l() const { return d_phi_l; };
  const std::vector<double> &get_phi_r() const { return d_phi_r; };

private:
  std::vector<double> d_phi_l;
  std::vector<double> d_phi_r;
};

class FETypeExteriorBdryNetwork {
public:
  explicit FETypeExteriorBdryNetwork()
    : d_phi(3) {
    d_phi[0] = 1;
    d_phi[1] = NAN;
    d_phi[2] = NAN;
  }

  /// Assume: (v)----e----
  void reinit(const Vertex & v, const Edge &e) {
    const double point_left = e.get_vertex_neighbors()[0] == v.get_id() ? -1 : +1;
    d_phi[1] = point_left;
    d_phi[2] = legendre2(point_left);
  };

  const std::vector<double> &get_phi() const { return d_phi; };

private:
  std::vector<double> d_phi;
};
*/

class FETypeNetwork {
public:
  explicit FETypeNetwork(QuadratureFormula qf)
      : d_qf(std::move(qf)),
        d_phi(1, std::vector<double>(d_qf.size())),
        d_dphi(1, std::vector<double>(d_qf.size())),
        d_JxW(d_qf.size()) {
    for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
      d_phi[0][qp] = 1;
      d_dphi[0][qp] = 0;
      d_JxW[qp] = NAN;
    }
  }

  void reinit(const Edge &e) {
    auto length = e.get_length();

    for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
      d_JxW[qp] = d_qf.ref_weights[qp] * length / 2.;
    }
  };

  const std::vector<std::vector<double>> &get_phi() const { return d_phi; };

  const std::vector<std::vector<double>> &get_dphi() const { return d_dphi; };

  const std::vector<double> &get_JxW() const { return d_JxW; };

private:
  QuadratureFormula d_qf;

  std::vector<std::vector<double>> d_phi;

  std::vector<std::vector<double>> d_dphi;

  std::vector<double> d_JxW;
};

class FETypeInnerBdryNetwork {
public:
  explicit FETypeInnerBdryNetwork()
    : d_phi_r(1), d_phi_l(1) {
    d_phi_l[0] = 1;
    d_phi_r[0] = 1;
  }

  /// Assume: ----e_left----(v)----e_right----
  void reinit(const Vertex & v, const Edge &e_left, const Edge &e_right) { };

  const std::vector<double> &get_phi_l() const { return d_phi_l; };
  const std::vector<double> &get_phi_r() const { return d_phi_r; };

private:
  std::vector<double> d_phi_l;
  std::vector<double> d_phi_r;
};

class FETypeExteriorBdryNetwork {
public:
  explicit FETypeExteriorBdryNetwork()
    : d_phi(1), d_n(1) {
    d_phi[0] = 1;
  }

  /// Assume: (v)----e----
  void reinit(const Vertex & v, const Edge &e) {
    d_n = (e.get_vertex_neighbors()[0] == v.get_id()) ? -1 : +1;
  };

  const std::vector<double> &get_phi() const { return d_phi; };
  const double &get_normal() const { return d_n; };

private:
  std::vector<double> d_phi;
  double d_n;
};

} // namespace macrocirculation

#endif //TUMORMODELS_FE_VALUES_H
