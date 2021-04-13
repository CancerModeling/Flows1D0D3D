////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_VESSEL_FORMULAS_HPP
#define TUMORMODELS_VESSEL_FORMULAS_HPP

#include <cmath>
#include <gmm.h>

namespace macrocirculation {

/*! @brief Saves the material constants on a vessel. */
struct VesselParameters {
  VesselParameters() = default;

  VesselParameters(double G0, double A0, double rho)
      : G0(G0), A0(A0), rho(rho) {}

  double G0{};
  double A0{};
  double rho{};
};

/*! @brief Calculates Eq. (2.8) from "Multi-scale modeling of flow and transport processes in arterial networks and tissue".
 *
 * @param h0 The thickness of the blood vessel.
 * @param E Young's modulus.
 * @param nu Poisson ratio.
 * @param A0 The vessel area at p=0.
 * @return
 */
inline double calculate_G0(double h0, double E, double nu, double A0) {
  return std::sqrt(M_PI) * h0 * E / ((1 - nu * nu) * std::sqrt(A0));
}

/*! @brief Calculates c0. */
inline double calculate_c0(double G0, double rho, double A0) {
  return std::sqrt(G0 / (2 * rho));
}

/*! @brief Calculates the flow Q from the characteristics. */
inline double calculate_Q(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return 0.5 * (w2 - w1) * std::pow((w1 + w2) / (8 * c0), 4) * A0;
}

/*! @brief Calculates the derivative of Q w.r.t. w1, the back propagating characteristic. */
inline double calculate_diff_Q_w1(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return -0.5 * std::pow((w1 + w2) / (8 * c0), 4) * A0 +
         2 * A0 * (w2 - w1) * std::pow(w2 + w1, 3) / std::pow(8 * c0, 4);
}

/*! @brief Calculates the derivative of Q w.r.t. w2, the forward propagating characteristic. */
inline double calculate_diff_Q_w2(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return +0.5 * std::pow((w1 + w2) / (8 * c0), 4) * A0 +
         2 * A0 * (w2 - w1) * std::pow(w2 + w1, 3) / std::pow(8 * c0, 4);
}

/*! @brief Calculates the total pressure from the characteristics. */
inline double calculate_p_from_w1w2(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return 0.5 * rho * std::pow(0.5 * (w2 - w1), 2) + G0 * (std::pow((w1 + w2) / (8 * c0), 2) - 1);
}

inline double calculate_p_from_QA(double Q, double A, double G0, double rho, double A0) {
  return 0.5 * rho * std::pow(Q / A, 2) + G0 * (std::sqrt(A / A0) - 1);
}

inline double calculate_static_p(double A, double G0, double A0) {
  return G0 * (std::sqrt(A / A0) - 1);
}

/*! @brief Derivative of the total pressure with respect to the backward characteristic. */
inline double calculate_diff_p_w1(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return -0.25 * rho * (w2 - w1) + 2 * G0 * (w1 + w2) / std::pow((8 * c0), 2);
}

/*! @brief Derivative of the total pressure with respect to the forward characteristic. */
inline double calculate_diff_p_w2(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return +0.25 * rho * (w2 - w1) + 2 * G0 * (w1 + w2) / std::pow((8 * c0), 2);
}

/*! @brief Evaluates the back propagating wave from Q and A. */
inline double calculate_W1_value(const double Q, const double A, const double G0, const double rho, const double A0) {
  return -Q / A + 4 * std::sqrt(G0 / (2 * rho)) * std::pow(A / A0, 1. / 4.);
}

/*! @brief Evaluates the forward propagating wave from Q and A. */
inline double calculate_W2_value(double Q, double A, double G0, double rho, double A0) {
  return +Q / A + 4 * std::sqrt(G0 / (2 * rho)) * std::pow(A / A0, 1. / 4.);
}

/*! @brief Convertes the characteristic variables to flow Q and area A. */
inline void convert_w1w2_to_QA(double w1, double w2, const VesselParameters &p, double &Q, double &A) {
  const double c0 = calculate_c0(p.G0, p.rho, p.A0);
  A = std::pow((w1 + w2) / (8 * c0), 4) * p.A0;
  Q = 0.5 * (w2 - w1) * A;
}

/*! @brief Solves the system of equation
 *           W1(Q_up, A_up) = W1
 *           W2(Q_up, A_up) = W2
 *         for (Q_up, Q_up) at some inner vertex between to vessel segments.
 */
inline void
solve_W12(double &Q_up, double &A_up, const double W1, const double W2, const double G0, const double rho, const double A0) {
  const double in = 1. / 8. * std::sqrt(2. * rho / G0) * (W2 + W1);
  A_up = A0 * std::pow(in, 4);
  Q_up = A_up / 2. * (W2 - W1);
}

/*! @brief Inflow boundary condition modeling a heart beat.
 *         These values are from T. Koeppls doctoral thesis, Eq. (2.124)
 *
 * @param t Current time.
 * @return The current flow rate.
 */
class heart_beat_inflow {
public:
  explicit heart_beat_inflow(double amplitude = 4.85, double t_period = 1.0, double t_systole = 0.3)
      : d_amplitude(amplitude),
        d_t_period(t_period),
        d_t_systole(t_systole) {}

  double operator()(double t) const {
    const double t_in_period = t - std::floor(t / d_t_period);
    if (t_in_period < d_t_systole) {
      return d_amplitude * std::sin(M_PI * t_in_period / d_t_systole);
    } else if (t_in_period <= d_t_period + 1e-14) {
      return 0.;
    } else {
      throw std::runtime_error("unreachable code");
    }
  }

private:
  double d_amplitude;
  double d_t_period;
  double d_t_systole;
};

/*! @brief Solves for the forward propagating characteristic W2, given the
 *         back propagating characteristic W1 and a prescribed flow Q_star,
 *         with a Newton iteration.
 *
 * @param W1      The backward propagating characteristic.
 * @param W2_init An initial guess for the forward propagating
 * @param Q_star  A prescribed flow.
 * @param A_0     The vessel area at zero pressure.
 * @param c_0     The velocity.
 * @return        The forward propagating wave w2.
 */
inline double solve_for_W2(const double W1, const double W2_init, const double Q_star, const double A_0, const double c_0) {
  double W2 = W2_init;
  double W2_prev = W2_init;
  const auto f = [=](double W2) { return (W2 - W1) / 2. * std::pow((W1 + W2) / (8 * c_0), 4) * A_0 - Q_star; };
  const auto f_prime = [=](double W2) {
    return (0.000610352 * A_0 * (W2 - 0.6 * W1) * pow((W2 + W1), 3)) / pow(c_0, 4);
  };
  for (std::size_t it = 0; it < 100; it += 1) {
    W2 += (-1) * f(W2) / f_prime(W2);
    if (std::abs(W2 - W2_prev) < 1e-16 || std::abs(W2 - W2_prev) < 1e-8 * W2_prev)
      break;
    W2_prev = W2;
  }
  return W2;
}

/*! @brief Solves for the backward propagating characteristic W1, given the
 *         forward propagating characteristic W2 and a prescribed flow Q_star,
 *         with a Newton iteration.
 *
 * @param W1_init An initial guess for the backward propagating characteristic.
 * @param W2      The forward
 * @param Q_star  A prescribed flow.
 * @param A_0     The vessel area at zero pressure.
 * @param c_0     The velocity.
 * @return        The forward propagating wave w2.
 */
inline double solve_for_W1(const double W1_init, const double W2, const double Q_star, const double A_0, const double c_0) {
  double W1 = W1_init;
  double W1_prev = W1_init;
  const auto f = [=](double W1) { return (W2 - W1) / 2. * std::pow((W1 + W2) / (8 * c_0), 4) * A_0 - Q_star; };
  const auto f_prime = [=](double W1) {
    return (-0.0001220703125 * A_0 * (5 * W1 - 3 * W2) * pow((W2 + W1), 3)) / pow(c_0, 4);
  };
  for (std::size_t it = 0; it < 100; it += 1) {
    W1 += (-1) * f(W1) / f_prime(W1);
    if (std::abs(W1 - W1_prev) < 1e-16 || std::abs(W1 - W1_prev) < 1e-8 * W1_prev)
      break;
    W1_prev = W1;
  }
  return W1;
}

/*! @brief Assembles the inflow boundary condition.
 *
 * @param Q  The value of Q inside the cell.
 * @param A  The value of A inside the cell.
 * @param in True, if the vessel points towards the vertex, false if it points away.
 * @param Q_star The Q value at the boundary
 * @param G0  TODO:
 * @param rho The blood density.
 * @param A0  The area at p=0.
 * @return
 */
inline double assemble_in_flow(const double Q,
                               const double A,
                               const bool in,
                               const double Q_star,
                               const double G0,
                               const double rho,
                               const double A0) {
  const double c0 = sqrt(G0 / (2 * rho));
  double W1 = calculate_W1_value(Q, A, G0, rho, A0);
  double W2 = calculate_W2_value(Q, A, G0, rho, A0);
  if (in)
    W1 = solve_for_W1(W1, W2, Q_star, A0, c0);
  else
    W2 = solve_for_W2(W1, W2, Q_star, A0, c0);
  const double A_star = A0 * std::pow(1. / (8 * c0) * (W2 + W1), 4);
  return A_star;
}

/*! @brief Evaluates the bifurcation equation at a vertex, for n vessels meeting at that point.
 *         For three vessels a, b and c, the bifurcation equation is given by
 *              B = (s_a*Q(a) + s_b*Q(b) + s_c*Q(c), p(a) - p(b), p(a) - p(c)),
 *         where s_a = -1 if the vessel a points towards the vertex and s_a = +1 otherwise.
 *
 * @param w1 The back propagating characteristic for vessel 1 to n.
 * @param w2 The forward propagating characteristic for vessel 1 to n.
 * @param p  The vessel parameters for vessel 1 to n.
 * @param in True if a vessel points towards the vertex, false if it points away.
 * @param out  A vector containing the three equations of the bifurcation system.
 */
inline void nfurcation_equation(const std::vector<double> &w1,
                                const std::vector<double> &w2,
                                const std::vector<VesselParameters> &p,
                                const std::vector<bool> &in,
                                std::vector<double> &out) {
  const size_t num_vessels = w1.size();

  assert(num_vessels == w2.size());
  assert(num_vessels == p.size());
  assert(num_vessels == in.size());
  assert(num_vessels == out.size());
  assert(num_vessels > 1);

  // flow equation (all flows add up to zero):
  double eq_Q = 0;
  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
    eq_Q += (in[vessel_idx] ? -1 : +1) * calculate_Q(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
  out[0] = eq_Q;

  // pressure equation (all pressures are equal):
  for (size_t vessel_idx = 1; vessel_idx < num_vessels; vessel_idx += 1)
    out[vessel_idx] = calculate_p_from_w1w2(w1[0], w2[0], p[0].G0, p[0].rho, p[0].A0) -
                      calculate_p_from_w1w2(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
}

/*! @brief Evaluates the jacobian nxn of the bifurcation equation system at a vertex,
 *         where n is the number of vessels meeting at one point.
 *         The evaluation depends on which vessels point towards the vertex.
 *         If a points towards the vertex and b and c point away, then we get
 *              [[ dQ(a)/dw_1, dQ(b)/dw_2, dQ(c)/dw_2],
 *               [ dp(a)/dw_1, dp(b)/dw_2, 0         ],
 *               [ dp(a)/dw_1, 0         , dp(c)/dw_2]]
 *         If the vessel point differently other derivatives are used.
 *
 * @param w1 The back propagating characteristic for vessel 1 to n.
 * @param w2 The forward propagating characteristic for vessel  1 to n.
 * @param p  The vessel parameters for vessel 1 to n.
 * @param in True if a vessel points towards the vertex, false if it points away.
 * @param J  The nxn Jacobian of the bifurcation equation system.
 */
inline void nfurcation_equation_jacobian(const std::vector<double> &w1,
                                         const std::vector<double> &w2,
                                         const std::vector<VesselParameters> &p,
                                         const std::vector<bool> &in,
                                         gmm::row_matrix<gmm::wsvector<double>> &J) {
  const size_t num_vessels = w1.size();

  assert(num_vessels == w2.size());
  assert(num_vessels == p.size());
  assert(num_vessels == in.size());
  assert(J.ncols() == num_vessels && J.nrows() == num_vessels);
  assert(num_vessels > 1);

  // line 0, Q derivative:
  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
    if (in[vessel_idx])
      J[0][vessel_idx] = -calculate_diff_Q_w1(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
    else
      J[0][vessel_idx] = +calculate_diff_Q_w2(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
  }

  // line 1-(n-1), pressure derivatives:
  for (size_t vessel_idx = 1; vessel_idx < num_vessels; vessel_idx += 1) {
    if (in[0])
      J[vessel_idx][0] = +calculate_diff_p_w1(w1[0], w2[0], p[0].G0, p[0].rho, p[0].A0);
    else
      J[vessel_idx][0] = +calculate_diff_p_w2(w1[0], w2[0], p[0].G0, p[0].rho, p[0].A0);

    if (in[vessel_idx])
      J[vessel_idx][vessel_idx] = -calculate_diff_p_w1(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
    else
      J[vessel_idx][vessel_idx] = -calculate_diff_p_w2(w1[vessel_idx], w2[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
  }
}

/*! @brief Takes the unknown characteristic quantities from a vector and
 *         copies them into the given values of the characteristics.
 *         E.g. if a vessel points towards a given vertex (in[vessel_idx] = true),
 *         then vec[vessel_idx] is copied to the unknown w1[vessel_idx], while w2[vessel_idx remains unchanged.
 *         This applies to all 0 <= vessel_idx < n.
 *
 * @param w1 The back propagating characteristic for vessel 1 to n.
 * @param w2 The forward propagating characteristic for vessel 1 to n.
 * @param in True if a vessel points towards the vertex, false if it points away.
 * @param vec  The vector of length n, whose quantities we want to transfer to w{1,2}.
 */
inline void extract_vector(std::vector<double> &w1,
                           std::vector<double> &w2,
                           const std::vector<bool> &in,
                           const std::vector<double> &vec) {
  const size_t num_vessels = w1.size();

  assert(num_vessels == w2.size());
  assert(num_vessels == in.size());
  assert(num_vessels == vec.size());
  assert(num_vessels > 1);

  // line 0, Q derivative:
  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
    if (in[vessel_idx])
      w1[vessel_idx] = vec[vessel_idx];
    else
      w2[vessel_idx] = vec[vessel_idx];
  }
}

/*! @brief Copies the unknown characteristic quantities into a vector.
 *         E.g. if vessel a points towards a given vertex (in[vessel_idx] = true),
 *         then w1[vessel_idx] is copied to vec[vessel_idx], while w2[vessel_idx] is not used at all.
 *         This applies to all 0 <= vessel_idx < num-vessels.
 *
 * @param w1 The back propagating characteristic for vessel 1 to n.
 * @param w2 The forward propagating characteristic for vessel 1 to n.
 * @param in True if vessel a points towards the vertex, false if it points away.
 * @param vec The vector, which we want to update with w{1,2}.
 */
inline void fill_vector(const std::vector<double> &w1,
                        const std::vector<double> &w2,
                        const std::vector<bool> &in,
                        std::vector<double> &vec) {
  const size_t num_vessels = w1.size();

  assert(num_vessels == w2.size());
  assert(num_vessels == in.size());
  assert(num_vessels == vec.size());
  assert(num_vessels > 1);

  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
    if (in[vessel_idx])
      vec[vessel_idx] = w1[vessel_idx];
    else
      vec[vessel_idx] = w2[vessel_idx];
  }
}

/*! @brief Solves the nfurcation equations at the intersection of n-vessels,
 *         and calculates the upwinded values for the respective vessels,
 *         taking into consideration the vessel orientation.
 *         Internally a Newton iteration is used to solve the system.
 *
 * @param w1   The back propagating characteristics for vessel 1 to n.
 * @param w2   The forward propagating characteristics for vessel 1 to n.
 * @param p    The vessel parameters for vessel 1 to n.
 * @param in   True if a vessel points towards the vertex, false if it points away.
 * @param Q_up The upwinded flow Q at vessel tips.
 * @param A_up The upwinded area A at vessel tips.
 */
inline std::size_t solve_at_nfurcation(const std::vector<double> &Q,
                                       const std::vector<double> &A,
                                       const std::vector<VesselParameters> &p,
                                       const std::vector<bool> &in,
                                       std::vector<double> &Q_up,
                                       std::vector<double> &A_up) {
  const size_t num_vessels = Q.size();

  assert(num_vessels == A.size());
  assert(num_vessels == p.size());
  assert(num_vessels == in.size());
  assert(num_vessels == Q_up.size());
  assert(num_vessels == A_up.size());
  assert(num_vessels > 1);

  const std::size_t max_iter = 1000;

  // calculate w1, w2 at all the bifurcations
  std::vector<double> w1(num_vessels, 0);
  std::vector<double> w2(num_vessels, 0);
  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
    w1[vessel_idx] = calculate_W1_value(Q[vessel_idx], A[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
    w2[vessel_idx] = calculate_W2_value(Q[vessel_idx], A[vessel_idx], p[vessel_idx].G0, p[vessel_idx].rho, p[vessel_idx].A0);
  }

  // solve for w1, w2 with a newton iteration
  gmm::row_matrix<gmm::wsvector<double>> Jac(num_vessels, num_vessels);
  std::vector<double> x(num_vessels, 0);
  std::vector<double> f(num_vessels, 0);
  std::vector<double> delta_x(num_vessels, 0);

  fill_vector(w1, w2, in, x);
  std::size_t it = 0;
  for (; it < max_iter; it += 1) {
    nfurcation_equation_jacobian(w1, w2, p, in, Jac);
    nfurcation_equation(w1, w2, p, in, f);
    gmm::lu_solve(Jac, delta_x, f);
    gmm::add(x, gmm::scaled(delta_x, -1), x);
    extract_vector(w1, w2, in, x);

    if (gmm::vect_norm2(delta_x) < 1e-14 || gmm::vect_norm2(delta_x) < 1e-14 * gmm::vect_norm2(x))
      break;
  }

  // copy back into the input variables
  for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
    convert_w1w2_to_QA(w1[vessel_idx], w2[vessel_idx], p[vessel_idx], Q_up[vessel_idx], A_up[vessel_idx]);

  return it;
}

} // namespace macrocirculation

#endif //TUMORMODELS_VESSEL_FORMULAS_HPP
