////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_FE_TYPE_HPP
#define TUMORMODELS_FE_TYPE_HPP

#include <utility>
#include <vector>

#include "graph_storage.hpp"

namespace macrocirculation {

// Collection of legendre polynomials
template<std::size_t degree>
constexpr double legendre(double);
template<>
constexpr double legendre<0>(double) {
  return 1.;
}
template<>
constexpr double legendre<1>(double x) {
  return x;
}
template<>
constexpr double legendre<2>(double x) {
  return 0.5 * (3 * x * x - 1);
}
template<>
constexpr double legendre<3>(double x) {
  return 0.5 * (5 * x * x * x - 3 * x);
}

// Collection of their derivatives
template<std::size_t degree>
constexpr double diff_legendre(double);
template<>
constexpr double diff_legendre<0>(double) {
  return 0;
}
template<>
constexpr double diff_legendre<1>(double /*x*/) {
  return 1;
}
template<>
constexpr double diff_legendre<2>(double x) {
  return 3 * x;
}
template<>
constexpr double diff_legendre<3>(double x) {
  return 0.5 * (15 * x * x - 3);
}

/*! @brief Representation of a quadrature formula for the network. */
struct QuadratureFormula {
  std::vector<double> ref_points;
  std::vector<double> ref_weights;

  std::size_t size() const { return ref_points.size(); }
};

/*! @brief Creates a gauss quadrature formula of order 7. */
inline QuadratureFormula create_gauss4() {
  QuadratureFormula qf;
  qf.ref_points = {-0.861136311594053, -0.339981043584856, +0.339981043584856, +0.861136311594053};
  qf.ref_weights = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
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

struct EdgeBoundaryValues {
  double left;
  double right;
};

/*! @brief Class for evaluating the legendre shape functions on our network edges. */
class FETypeNetwork {
public:
  FETypeNetwork(QuadratureFormula qf, std::size_t degree);

  /*! @brief Updates the shape function values on the given edge. */
  void reinit(double length);

  /*! @returns Returns a list of basis functions evaluated at the quadrature points.
   *           The list is structured by phi[ <shape-function-index> ][ <quadrature-point-index> ],
   *           which is the same as in libmesh.
   */
  const std::vector<std::vector<double>> &get_phi() const { return d_phi; };

  /*! @returns Returns a list of basis functions evaluated at the boundary poitns.
   *           The list is structured by phi[ <shape-function-index> ][ <boundary-index> ],
   *           where the left boundary has index 0 and the right boundary an index of 1.
   */
  const std::vector<std::vector<double>> &get_phi_boundary() const { return d_phi_boundary; };

  /*! @returns Returns a list of derivative of basis functions evaluated at the quadrature points.
   *           The list is structured by dphi[ <shape-function-index> ][ <quadrature-point-index> ],
   *           which is the same as in libmesh.
   */
  const std::vector<std::vector<double>> &get_dphi() const { return d_dphi; };

  /*! @returns Returns a list of quadrature weights at each quadrature point. */
  const std::vector<double> &get_JxW() const { return d_JxW; };

  /*! @brief Evaluates the function with the given dof values at the quadratures points. */
  void evaluate_dof_at_quadrature_points(const std::vector<double> &dof_values,
                                         std::vector<double> &quadrature_point_values) const;

  /*! @brief Evaluates the finite element basis functions on the interval \f$[-1,+1]\f$ at \f$s\in[-1,+1]\f$. */
  static double evaluate_dof(const std::vector<double> &dof_values, double s);

  /*! @brief Evaluates the function with the given dof values at the quadratures points. */
  EdgeBoundaryValues evaluate_dof_at_boundary_points(const std::vector<double> &dof_values) const;

  /*! @brief Returns the dof values for the linear function with value v0 at the left and v1 at the right. */
  void interpolate_linear_function(double v0, double v1, std::vector<double> &dof_values) const;

  /*! @returns The number of quadrature point for which this finite-element was constructed. */
  std::size_t num_quad_points() const;

  /*! @returns The degree of the legendre basis. */
  std::size_t get_degree() const;

  /*! @returns The quadrature formula for the given finite element type. */
  QuadratureFormula get_quadrature_formula() const;

private:
  QuadratureFormula d_qf;

  std::size_t d_degree;

  std::vector<std::vector<double>> d_phi;

  std::vector<std::vector<double>> d_phi_boundary;

  std::vector<std::vector<double>> d_dphi;

  std::vector<double> d_JxW;
};

class QuadraturePointMapper {
public:
  explicit QuadraturePointMapper(const QuadratureFormula &qf);

  void reinit(double s_left, double s_right);

  const std::vector<double> &get_quadrature_points() const { return d_quadrature_points; };

private:
  QuadratureFormula d_qf;

  std::vector<double> d_quadrature_points;
};

} // namespace macrocirculation

#endif