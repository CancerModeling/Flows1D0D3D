////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#include "stochastic_noise_assembly.hpp"

namespace util {

StochasticNoiseAssembly::StochasticNoiseAssembly(unsigned int num_eigenfunctions,
                                                 unsigned int seed,
                                                 double scale,
                                                 double length,
                                                 double lower_bound,
                                                 double upper_bound)
    : d_num_terms_bound(num_eigenfunctions),
      d_generator(seed),
      d_scale(scale),
      d_length(length),
      d_lower_bound(lower_bound),
      d_upper_bound(upper_bound),
      d_avg(0.),
      d_cached_rhs(nullptr),
      d_reassemble(true) {}

void StochasticNoiseAssembly::assemble(BaseAssembly &assembly, BaseAssembly &total_assembly) const {
  // create cached vector if not present yet and reassemble vector only if necessary
  if (d_cached_rhs == nullptr || d_reassemble) {
    d_cached_rhs = assembly.d_sys.rhs->zero_clone();
    assemble(assembly, total_assembly, *d_cached_rhs);
    d_reassemble = false;
  }

  // add stochastic contributions
  assembly.d_sys.rhs->add(*d_cached_rhs);
}

void StochasticNoiseAssembly::assemble(BaseAssembly &assembly, BaseAssembly &total_assembly, libMesh::NumericVector<Number> &rhs) const {
  // if the stochastic scaling factor is too small, we skip assembly,
  // since evaluation of the noise at each quadrature point is not really cheap.
  if (std::abs(d_scale) < 1e-14 || d_num_terms_bound == 0)
    return;

  const auto &quad_points = assembly.d_fe->get_xyz();
  const auto &phi = assembly.d_phi;
  const auto &JxW = assembly.d_JxW;
  const auto &dof_indices_sys = assembly.d_dof_indices_sys_var[0];

  libMesh::DenseVector<double> f;

  for (const auto &elem : assembly.d_mesh.active_local_element_ptr_range()) {
    assembly.init_dof(elem);
    assembly.init_fe(elem);
    total_assembly.init_dof(elem);
    for (unsigned int qp = 0; qp < assembly.d_qrule.n_points(); qp++) {
      Real field_value_old = 0;
      Real total_field_value_old = 0;
      for (unsigned int i = 0; i < phi.size(); i++) {
        field_value_old += phi[i][qp] * assembly.get_old_sol_var(i, 0);
        total_field_value_old += phi[i][qp] * total_assembly.get_old_sol_var(i, 0);
      }
      const auto value_at_qp = eval_eigenfunctions_at_quadrature_point(quad_points[qp], field_value_old, total_field_value_old);
      f.resize(phi.size());
      for (unsigned int i = 0; i < phi.size(); i++) {
        f(i) = JxW[qp] * d_scale * value_at_qp * phi[i][qp];
      }
      rhs.add_vector(f, dof_indices_sys);
    }
  }
}

void StochasticNoiseAssembly::calculate_new_stochastic_coefficients(const double dt) {

  // return if we do not want stochastic term
  if (std::abs(d_scale) < 1e-14 || d_num_terms_bound == 0)
    return;

  std::normal_distribution<double> gaussian(0, std::sqrt(dt));

  d_stochastic_coefficients.resize(num_basis_functions());

  for (int i = 0; i < num_basis_functions(); i += 1)
    d_stochastic_coefficients[i] = gaussian(d_generator);

  // set the flag to force a reassembly
  d_reassemble = true;
}

void StochasticNoiseAssembly::calculate_new_stochastic_coefficients(const double dt, BaseAssembly &assembly, BaseAssembly &total_assembly) {

  // return if we do not want stochastic term
  if (std::abs(d_scale) < 1e-14 || d_num_terms_bound == 0)
    return;

  std::normal_distribution<double> gaussian(0, std::sqrt(dt));

  d_stochastic_coefficients.resize(num_basis_functions());

  for (int i = 0; i < num_basis_functions(); i += 1)
    d_stochastic_coefficients[i] = gaussian(d_generator);

  // set the flag to force a reassembly
  d_reassemble = true;

  // compute average
  double s = 0.;
  double h = 0.;
  d_avg = 0.; // reset to zero so that we get purely H \Delta W term in below

  // if the stochastic scaling factor is too small, we skip assembly,
  // since evaluation of the noise at each quadrature point is not really cheap.
  if (std::abs(d_scale) < 1e-14 || d_num_terms_bound == 0)
    return;

  const auto &quad_points = assembly.d_fe->get_xyz();
  const auto &phi = assembly.d_phi;
  const auto &JxW = assembly.d_JxW;
  const auto &dof_indices_sys = assembly.d_dof_indices_sys_var[0];

  for (const auto &elem : assembly.d_mesh.active_local_element_ptr_range()) {
    assembly.init_dof(elem);
    assembly.init_fe(elem);
    total_assembly.init_dof(elem);
    for (unsigned int qp = 0; qp < assembly.d_qrule.n_points(); qp++) {
      Real field_value_old = 0;
      Real total_field_value_old = 0;
      for (unsigned int i = 0; i < phi.size(); i++) {
        field_value_old += phi[i][qp] * assembly.get_old_sol_var(i, 0);
        total_field_value_old += phi[i][qp] * total_assembly.get_old_sol_var(i, 0);
      }

      // heaviside function
      Real heavyside_value = util::heaviside(field_value_old - d_lower_bound) * util::heaviside(d_upper_bound - field_value_old) *
                             util::heaviside(total_field_value_old - d_lower_bound) * util::heaviside(d_upper_bound - total_field_value_old);

      // get H \Delta W at quad point
      const auto value_at_qp = eval_eigenfunctions_at_quadrature_point(quad_points[qp], field_value_old, total_field_value_old);

      // compute average of H \Delta W and H
      for (unsigned int i = 0; i < phi.size(); i++) {
        s += JxW[qp] * d_scale * value_at_qp * phi[i][qp];
        h += JxW[qp] * heavyside_value * phi[i][qp];
      }
    }
  }

  // communicate with other processors
  double s_global = 0.;
  double h_global = 0.;
  MPI_Allreduce(&s, &s_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&h, &h_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // set average
  if (h_global > 1.e-10)
    d_avg = s_global / h_global;
}

double StochasticNoiseAssembly::eval_eigenfunctions_at_quadrature_point(const Point &p, const double field_value, const double total_field_value) const {

  if (d_stochastic_coefficients.size() != num_basis_functions())
    throw std::runtime_error("not enough stochastic coefficients precalculated");

  auto heavyside_value = util::heaviside(field_value - d_lower_bound) * util::heaviside(d_upper_bound - field_value) *
                               util::heaviside(total_field_value - d_lower_bound) * util::heaviside(d_upper_bound - total_field_value);

  // quit if we are not in the interval
  if (heavyside_value < 1e-12)
    return 0;

  if (total_field_value <= d_lower_bound || total_field_value >= d_upper_bound)
    return 0;

  double acc = 0;

  const double normalization = std::sqrt(8 / std::pow(d_length, 3));
  std::size_t idx = 0;
  for (int kx = 0; kx < d_num_terms_bound; kx += 1) {
    const double mul_x = std::cos(kx * M_PI / d_length * p(0));
    for (int ky = 0; ky < d_num_terms_bound-kx; ky += 1) {
      const double mul_y = std::cos(ky * M_PI / d_length * p(1));
      for (int kz = 0; kz < d_num_terms_bound-kx-ky; kz += 1) {
        const double mul_z = std::cos(kz * M_PI / d_length * p(2));
        acc += d_stochastic_coefficients[idx++] * heavyside_value * normalization * mul_x * mul_y * mul_z;
      }
    }
  }

  // subtract average and return (d_avg will be zero if user has set d_hyp_subtract_avg_stoch to false in input)
  return acc - heavyside_value * d_avg;
}

std::size_t StochasticNoiseAssembly::num_basis_functions() const {
  return d_num_terms_bound*(d_num_terms_bound+1)*(d_num_terms_bound+2)/6.;
}

} // namespace util
