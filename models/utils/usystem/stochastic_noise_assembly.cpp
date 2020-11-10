#include "stochastic_noise_assembly.hpp"

namespace util {

StochasticNoiseAssembly::StochasticNoiseAssembly(unsigned int num_eigenfunctions,
                                                 unsigned int seed,
                                                 double scale,
                                                 double length,
                                                 double lower_bound,
                                                 double upper_bound)
    : d_num_eigenfunctions(num_eigenfunctions),
      d_generator(seed),
      d_scale(scale),
      d_length(length),
      d_lower_bound(lower_bound),
      d_upper_bound(upper_bound) {}

void StochasticNoiseAssembly::assemble(BaseAssembly &assembly) const {
  const auto &quad_points = assembly.d_fe->get_xyz();
  const auto &phi = assembly.d_phi;
  const auto &JxW = assembly.d_JxW;
  const auto &dof_indices_sys = assembly.d_dof_indices_sys;

  libMesh::DenseVector<double> f;

  for (const auto &elem : assembly.d_mesh.active_local_element_ptr_range()) {
    assembly.init_dof(elem);
    for (unsigned int qp = 0; qp < assembly.d_qrule.n_points(); qp++) {
      Real field_value_old = 0;
      for (unsigned int i = 0; i < phi.size(); i++) {
        field_value_old += phi[i][qp] * assembly.get_old_sol_var(i, 0);
      }
      const auto value_at_qp = eval_eigenfunctions_at_quadrature_point(quad_points[qp], field_value_old);
      f.resize(phi.size());
      for (unsigned int i = 0; i < phi.size(); i++) {
        f(i) = JxW[qp] * d_scale * value_at_qp * phi[i][qp];
      }
      assembly.d_sys.rhs->add_vector(f, dof_indices_sys);
    }
  }
}

void StochasticNoiseAssembly::calculate_new_stochastic_coefficients(const double dt) {
  const int N = std::ceil(std::cbrt(d_num_eigenfunctions));

  std::normal_distribution<double> gaussian(0, std::sqrt(dt));

  if (std::abs(static_cast<double>(N * N * N) - d_num_eigenfunctions) > 1e-14)
    throw std::runtime_error("the number of eigenfunctions must be of the form N^3");

  d_stochastic_coefficients.resize(N * N * N);

  for (int i = 0; i < N * N * N; i += 1)
    d_stochastic_coefficients[i] = gaussian(d_generator);
}

double StochasticNoiseAssembly::eval_eigenfunctions_at_quadrature_point(const Point &p, const double field_value) const {
  const int N = std::ceil(std::cbrt(d_num_eigenfunctions));

  if (std::abs(static_cast<double>(N * N * N) - d_num_eigenfunctions) > 1e-14)
    throw std::runtime_error("the number of eigenfunctions must be of the form N^3");

  if (d_stochastic_coefficients.size() != N * N * N)
    throw std::runtime_error("not enough stochastic coefficients precalculated");

  // quit if we are not in the interval
  if (field_value <= d_lower_bound || field_value >= d_upper_bound)
    return 0;

  double acc = 0;

  const double normalization = std::sqrt(8 / std::pow(d_length, 3));
  for (int kx = 0; kx < N; kx += 1) {
    const double mul_x = std::cos(kx * M_PI / d_length * p(0));
    for (int ky = 0; ky < N; ky += 1) {
      const double mul_y = std::cos(ky * M_PI / d_length * p(1));
      for (int kz = 0; kz < N; kz += 1) {
        const double mul_z = std::cos(kz * M_PI / d_length * p(2));
        acc += d_stochastic_coefficients[N * N * kx + N * ky + kz] * normalization * mul_x * mul_y * mul_z;
      }
    }
  }
  return acc;
}

} // namespace util
