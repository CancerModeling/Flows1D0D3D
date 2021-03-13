////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "fe_type.hpp"

#include <cassert>
#include <cmath>

namespace macrocirculation {

FETypeNetwork::FETypeNetwork(QuadratureFormula qf, std::size_t degree)
    : d_qf(std::move(qf)),
      d_degree(degree),
      d_phi(d_degree + 1, std::vector<double>(d_qf.size())),
      d_phi_boundary(d_degree + 1, std::vector<double>(2)),
      d_dphi(d_degree + 1, std::vector<double>(d_qf.size())), d_JxW(d_qf.size()) {
  // we only have lagrange polynomials up to order three
  assert(d_degree < 4);

  // phi inside of cell
  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    d_phi[0][qp] = 1;
    d_dphi[0][qp] = 0;

    if (d_degree > 0) {
      d_phi[1][qp] = legendre<1>(d_qf.ref_points[qp]);
      d_dphi[1][qp] = NAN;
    }

    if (d_degree > 1) {
      d_phi[2][qp] = legendre<2>(d_qf.ref_points[qp]);
      d_dphi[2][qp] = NAN;
    }

    if (d_degree > 2) {
      d_phi[3][qp] = legendre<3>(d_qf.ref_points[qp]);
      d_dphi[3][qp] = NAN;
    }

    d_JxW[qp] = NAN;
  }

  // phi on boundary
  {
    d_phi_boundary[0][0] = 1;
    d_phi_boundary[0][1] = 1;

    if (d_degree > 0) {
      d_phi_boundary[1][0] = legendre<1>(-1.);
      d_phi_boundary[1][1] = legendre<1>(+1.);
    }

    if (d_degree > 1) {
      d_phi_boundary[2][0] = legendre<2>(-1.);
      d_phi_boundary[2][1] = legendre<2>(+1.);
    }

    if (d_degree > 2) {
      d_phi_boundary[3][0] = legendre<3>(-1.);
      d_phi_boundary[3][1] = legendre<3>(+1.);
    }
  }
}

void FETypeNetwork::reinit(double length) {
  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    if (d_degree > 0)
      d_dphi[1][qp] = diff_legendre<1>(d_qf.ref_points[qp]) * 2. / length;

    if (d_degree > 1)
      d_dphi[2][qp] = diff_legendre<2>(d_qf.ref_points[qp]) * 2. / length;

    if (d_degree > 2)
      d_dphi[3][qp] = diff_legendre<3>(d_qf.ref_points[qp]) * 2. / length;

    d_JxW[qp] = d_qf.ref_weights[qp] * length / 2.;
  }
}

void FETypeNetwork::evaluate_dof_at_quadrature_points(const std::vector<double> &dof_values,
                                                      std::vector<double> &quadrature_point_values) const {
  assert(dof_values.size() == d_degree + 1);
  assert(quadrature_point_values.size() == d_qf.size());

  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    quadrature_point_values[qp] = 0;
    for (std::size_t i = 0; i < d_degree + 1; i += 1) {
      quadrature_point_values[qp] += d_phi[i][qp] * dof_values[i];
    }
  }
}

std::size_t FETypeNetwork::num_quad_points() const {
  return d_qf.size();
}

std::size_t FETypeNetwork::get_degree() const {
  return d_degree;
}

EdgeBoundaryValues FETypeNetwork::evaluate_dof_at_boundary_points(const std::vector<double> &dof_values) const {
  EdgeBoundaryValues values{0, 0};

  values.left += legendre<0>(-1) * dof_values[0];
  values.right += legendre<0>(+1) * dof_values[0];

  if (d_degree > 0) {
    values.left += legendre<1>(-1) * dof_values[1];
    values.right += legendre<1>(+1) * dof_values[1];
  }

  if (d_degree > 1) {
    values.left += legendre<2>(-1) * dof_values[2];
    values.right += legendre<2>(+1) * dof_values[2];
  }

  if (d_degree > 2) {
    values.left += legendre<3>(-1) * dof_values[3];
    values.right += legendre<3>(+1) * dof_values[3];
  }

  return values;
}

QuadraturePointMapper::QuadraturePointMapper(const QuadratureFormula &qf)
    : d_qf(qf), d_quadrature_points(qf.size()) {}

void QuadraturePointMapper::reinit(double s_left, double s_right) {
  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    const double r = 0.5 * (d_qf.ref_points[qp] + 1);
    d_quadrature_points[qp] = s_left * (1 - r) + s_right * r;
  }
}

} // namespace macrocirculation
