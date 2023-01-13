////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_UPWINDING_FORMULAS_LINEARIZED_FLOW_HPP
#define TUMORMODELS_UPWINDING_FORMULAS_LINEARIZED_FLOW_HPP

#include <Eigen/Dense>
#include <cassert>

#include "vessel_formulas.hpp"

namespace macrocirculation {

namespace linearized {

/*! @brief Calculates the upwinding at an inner boundary.
 *
 * @param alpha The \f$ \sqrt{\frac{C}{L}} \f$ factor for the linearized characteristics.
 * @param p_l   The pressure at the boundary of the left cell.
 * @param q_l   The flux at the boundary of the left cell.
 * @param p_r   The pressure at the boundary of the right cell.
 * @param q_r   The flux at the boundary of the right cell.
 * @param p_up  The upwinded pressure.
 * @param q_up  The upwinded flux.
 */
inline void inner_boundary(double alpha, double p_l, double q_l, double p_r, double q_r, double &p_up, double &q_up) {
  const double w1 = 0.5 * (-alpha * p_r  + q_r );
  const double w2 = 0.5 * (+alpha * p_l  + q_l );
  p_up = (w2 - w1) / alpha;
  q_up = w1 + w2;
}

/*! @brief Calculates the upwinding at an nfurcation.
 *
 * @param p         The boundary pressure values of each vessel at a common vertex.
 * @param q         The boundary flux values of each vessel at a common vertex.
 * @param params    The physical edge parameters near the nfurcation.
 * @param sigma     The vessel normals (+1 if the vessel points to the vertex, -1 if not).
 * @param p_up      The upwinded pressure values.
 * @param q_up      The upwinded flux values.
 */
inline void nfurcation_boundary(const std::vector<double> &p,
                                const std::vector<double> &q,
                                const std::vector<VesselParameters> &params,
                                const std::vector<double> &sigma,
                                std::vector<double> &p_up,
                                std::vector<double> &q_up) {
  // all vectors need to have the same size:
  assert(p.size() == q.size());
  assert(p_up.size() == p.size());
  assert(q_up.size() == p.size());
  assert(params.size() == p.size());
  assert(sigma.size() == p.size());

  const size_t N = p.size();

  // dof-ordering: (p_1, q_1, p_2, q_2, ... p_N, q_N)
  Eigen::MatrixXd mat(p.size() + q.size(), p.size() + q.size());
  Eigen::VectorXd rhs(p.size() + q.size());

  mat.setZero();

  // vector containing the \f$ \sqrt{ \frac{C}{L} } \f$ factors:
  std::vector<double> alpha;
  for (auto &param : params) {
    alpha.push_back(std::sqrt(linear::get_C(param) / linear::get_L(param)));
  }

  // constrain the upwinded fluxes to zero:
  for (int k = 0; k < static_cast<int>(N); k += 1) {
    mat(0, 2 * k + 1) = sigma[k];
  }
  rhs(0) = 0;

  // the characteristics should be equal
  for (int k = 0; k < static_cast<int>(N); k += 1) {
    mat(1 + k, 2 * k) = 0.5 * alpha[k] * sigma[k];
    mat(1 + k, 2 * k + 1) = 0.5;
    rhs(1 + k) = 0.5 * alpha[k] * sigma[k] * p[k] + 0.5 * q[k];
  }

  // the pressures should be equal
  for (int k = 1; k < static_cast<int>(N); k += 1) {
    mat(N + k, 2 * k) = 1.;
    mat(N + k, 2 * (k - 1)) = -1.;
    rhs(N + k) = 0.;
  }

  Eigen::VectorXd result = mat.fullPivLu().solve(rhs);

  for (int k = 0; k < static_cast<int>(N); k += 1) {
    p_up[k] = result[2 * k];
    q_up[k] = result[2 * k + 1];
  }
}

} // namespace linearized

} // namespace macrocirculation

#endif //TUMORMODELS_UPWINDING_FORMULAS_LINEARIZED_FLOW_HPP
