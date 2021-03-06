////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "fe_type_network.hpp"

#include <cassert>

namespace macrocirculation {

template<std::size_t DEGREE>
FETypeNetwork<DEGREE>::FETypeNetwork(QuadratureFormula qf)
    : d_qf(std::move(qf)),
      d_phi(DEGREE + 1, std::vector<double>(d_qf.size())),
      d_dphi(DEGREE + 1, std::vector<double>(d_qf.size())),
      d_JxW(d_qf.size()) {
  // we only have lagrange polynomials up to order three
  static_assert(DEGREE < 4);

  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    d_phi[0][qp] = 1;
    d_dphi[0][qp] = 0;

    if (DEGREE > 0) {
      d_phi[1][qp] = legendre1(d_qf.ref_points[qp]);
      d_dphi[1][qp] = NAN;
    }

    if (DEGREE > 1) {
      d_phi[2][qp] = legendre2(d_qf.ref_points[qp]);
      d_dphi[2][qp] = NAN;
    }

    if (DEGREE > 2) {
      d_phi[3][qp] = legendre3(d_qf.ref_points[qp]);
      d_dphi[3][qp] = NAN;
    }

    d_JxW[qp] = NAN;
  }
}

template<std::size_t DEGREE>
void FETypeNetwork<DEGREE>::reinit(const Edge &e) {
  auto length = e.get_length();

  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    if (DEGREE > 0)
      d_dphi[1][qp] = diff_legendre1(d_qf.ref_points[qp]) * 2. / length;

    if (DEGREE > 1)
      d_dphi[2][qp] = diff_legendre2(d_qf.ref_points[qp]) * 2. / length;

    if (DEGREE > 2)
      d_dphi[3][qp] = diff_legendre3(d_qf.ref_points[qp]) * 2. / length;

    d_JxW[qp] = d_qf.ref_weights[qp] * length / 2.;
  }
}

template<std::size_t DEGREE>
void FETypeNetwork<DEGREE>::evaluate_dof_at_quadrature_points(const std::vector<double> &dof_values, std::vector<double> &quadrature_point_values) const {
  assert(dof_values.size() == DEGREE + 1);
  assert(quadrature_point_values.size() == d_qf.size());

  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    quadrature_point_values[qp] = 0;
    for (std::size_t i = 0; i < DEGREE + 1; i += 1) {
      quadrature_point_values[qp] += d_phi[i][qp] * dof_values[i];
    }
  }
}

template<std::size_t DEGREE>
std::size_t FETypeNetwork<DEGREE>::num_quad_points() const {
  return d_qf.size();
}

template<std::size_t DEGREE>
FETypeInnerBdryNetwork<DEGREE>::FETypeInnerBdryNetwork()
    : d_phi_r(DEGREE + 1), d_phi_l(DEGREE + 1) {
  // we only have lagrange polynomials up to order three
  static_assert(DEGREE < 4);

  d_phi_l[0] = 1;
  d_phi_r[0] = 1;

  const double point_left = +1;
  const double point_right = -1;

  if (DEGREE > 0) {
    d_phi_l[1] = legendre1(point_left);
    d_phi_r[1] = legendre1(point_right);
  }

  if (DEGREE > 1) {
    d_phi_l[2] = legendre2(point_left);
    d_phi_r[2] = legendre2(point_right);
  }

  if (DEGREE > 2) {
    d_phi_l[3] = legendre3(point_left);
    d_phi_r[3] = legendre3(point_right);
  }
}

template<std::size_t DEGREE>
void FETypeInnerBdryNetwork<DEGREE>::reinit(const Vertex &v, const Edge &e_left, const Edge &e_right) {
  assert(e_left.get_vertex_neighbors()[1] == v.get_id());
  assert(e_right.get_vertex_neighbors()[0] == v.get_id());
};

template<std::size_t DEGREE>
EdgeBoundaryValues FETypeInnerBdryNetwork<DEGREE>::evaluate_dof_at_boundary_points(const std::vector<double> &dof_values) const {
  EdgeBoundaryValues values{0, 0};
  for (std::size_t i = 0; i < DEGREE + 1; i += 1) {
    values.left += d_phi_l[i] * dof_values[i];
    values.right += d_phi_r[i] * dof_values[i];
  }
  return values;
}

template<std::size_t DEGREE>
FETypeExteriorBdryNetwork<DEGREE>::FETypeExteriorBdryNetwork()
    : d_phi(DEGREE + 1), d_n(1) {
  // we only have lagrange polynomials up to order three
  static_assert(DEGREE < 4);

  d_phi[0] = 1;
  if (DEGREE > 0)
    d_phi[1] = NAN;
  if (DEGREE > 1)
    d_phi[2] = NAN;
  if (DEGREE > 3)
    d_phi[3] = NAN;
}

template<std::size_t DEGREE>
void FETypeExteriorBdryNetwork<DEGREE>::reinit(const Vertex &v, const Edge &e) {
  const double point = e.get_vertex_neighbors()[0] == v.get_id() ? -1 : +1;
  if (DEGREE > 0)
    d_phi[1] = legendre1(point);
  if (DEGREE > 1)
    d_phi[2] = legendre2(point);
  if (DEGREE > 2)
    d_phi[3] = legendre3(point);
  d_n = (e.get_vertex_neighbors()[0] == v.get_id()) ? -1 : +1;
};

template<std::size_t DEGREE>
double FETypeExteriorBdryNetwork<DEGREE>::evaluate_dof_at_boundary_points(const std::vector<double> &dof_values) const {
  double value = 0;
  for (std::size_t i = 0; i < DEGREE + 1; i += 1) {
    value += d_phi[i] * dof_values[i];
  }
  return value;
}

QuadraturePointMapper::QuadraturePointMapper(QuadratureFormula qf)
    : d_qf(qf),
      d_quadrature_points(qf.size()) {}

void QuadraturePointMapper::reinit(const Edge &e) {
  for (std::size_t qp = 0; qp < d_qf.size(); qp += 1) {
    const double s = 0.5 * (d_qf.ref_points[qp] + 1);
    d_quadrature_points[qp] = e.get_coordinate_v0() * (1 - s) + e.get_coordinate_v1() * s;
  }
}

// instantiations of the available templates to avoid lengthy recompiles:
// DG0:
template class FETypeNetwork<0>;
template class FETypeInnerBdryNetwork<0>;
template class FETypeExteriorBdryNetwork<0>;

// DG1:
template class FETypeNetwork<1>;
template class FETypeInnerBdryNetwork<1>;
template class FETypeExteriorBdryNetwork<1>;

// DG2:
template class FETypeNetwork<2>;
template class FETypeInnerBdryNetwork<2>;
template class FETypeExteriorBdryNetwork<2>;

// DG3:
template class FETypeNetwork<3>;
template class FETypeInnerBdryNetwork<3>;
template class FETypeExteriorBdryNetwork<3>;

} // namespace macrocirculation
