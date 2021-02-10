////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_FE_TYPE_NETWORK_H
#define TUMORMODELS_FE_TYPE_NETWORK_H

#include <utility>

#include "graph_storage.hpp"
#include "libmesh/libmesh.h"

namespace macrocirculation {

// Collection of legendre polynomials
constexpr double legendre1(double x) { return x; }
constexpr double legendre2(double x) { return 0.5 * (3 * x * x - 1); }
constexpr double legendre3(double x) { return 0.5 * (5 * x * x * x - 3 * x); }

// Collection of their derivatives
constexpr double diff_legendre1(double x) { return 1; }
constexpr double diff_legendre2(double x) { return 3 * x; }
constexpr double diff_legendre3(double x) { return 0.5 * (15 * x * x - 3); }

/*! @brief Representation of a quadrature formula for the network. */
struct QuadratureFormula {
  std::vector<double> ref_points;
  std::vector<double> ref_weights;

  std::size_t size() const { return ref_points.size(); }
};

/*! @brief Creates a gauss quadrature formula of order 7. */
inline QuadratureFormula create_gauss4() {
  QuadratureFormula qf;
  qf.ref_points = {
    -0.861136311594053,
    -0.339981043584856,
    +0.339981043584856,
    +0.861136311594053};
  qf.ref_weights = {
    0.347854845137454,
    0.652145154862546,
    0.652145154862546,
    0.347854845137454};
  return qf;
}

/*! @brief Creates a gauss quadrature formula of order 5. */
inline QuadratureFormula create_gauss3() {
  QuadratureFormula qf;
  qf.ref_points = {-0.774596669241, 0, 0.774596669241};
  qf.ref_weights = {+5. / 9., 8. / 9., +5. / 9.};
  return qf;
}

/*! @brief Creates tabulates the weights and points of the trapezoidal rule. */
inline QuadratureFormula create_trapezoidal_rule() {
  QuadratureFormula qf;
  qf.ref_points = {-1, +1};
  qf.ref_weights = {1, 1};
  return qf;
}

/*! @brief Creates tabulates the weights and points of the midpoint rule. */
inline QuadratureFormula create_midpoint_rule() {
  QuadratureFormula qf;
  qf.ref_points = {0};
  qf.ref_weights = {2.};
  return qf;
}

/*! @brief Class for evaluating the legendre shape functions on our network edges. */
template<std::size_t DEGREE>
class FETypeNetwork {
public:
  explicit FETypeNetwork(QuadratureFormula qf);

  /*! @brief Updates the shape function values on the given edge. */
  void reinit(const Edge &e);

  const std::vector<std::vector<double>> &get_phi() const { return d_phi; };

  const std::vector<std::vector<double>> &get_dphi() const { return d_dphi; };

  const std::vector<double> &get_JxW() const { return d_JxW; };

  /*! @brief Evaluates the function with the given dof values at the quadratures points. */
  void evaluate_dof_at_quadrature_points(const std::vector<double> &dof_values, std::vector<double> &quadrature_point_values) const;

  std::size_t num_quad_points() const;

private:
  QuadratureFormula d_qf;

  std::vector<std::vector<double>> d_phi;

  std::vector<std::vector<double>> d_dphi;

  std::vector<double> d_JxW;
};

/*! @brief Class for evaluating the legendre shape functions on our _inner_ network vertices. */
template<std::size_t DEGREE>
class FETypeInnerBdryNetwork {
public:
  FETypeInnerBdryNetwork();

  /*! @brief Updates the shape function values at the boundary.
   *         This class assumes that the edge e_left points towards the vertex v,
   *         while the edge e_right points away from the vertex.
   *         If this is not the case an exception is thrown.
   *
   * Geometric assumption: ----e_left----(v)----e_right----
   */
  void reinit(const Vertex &v, const Edge &e_left, const Edge &e_right);

  const std::vector<double> &get_phi_l() const { return d_phi_l; };
  const std::vector<double> &get_phi_r() const { return d_phi_r; };

private:
  std::vector<double> d_phi_l;
  std::vector<double> d_phi_r;
};

/*! @brief Class for evaluating the legendre shape functions on our exterior network vertices. */
template<std::size_t DEGREE>
class FETypeExteriorBdryNetwork {
public:
  FETypeExteriorBdryNetwork();

  /*! @brief Updates the shape function values on the given vertex. */
  void reinit(const Vertex &v, const Edge &e);

  const std::vector<double> &get_phi() const { return d_phi; };
  const double &get_normal() const { return d_n; };

  double evaluate_dof_at_boundary_points(const std::vector<double> &dof_values) const;

private:
  std::vector<double> d_phi;

  double d_n;
};

} // namespace macrocirculation

#endif //TUMORMODELS_FE_VALUES_H
