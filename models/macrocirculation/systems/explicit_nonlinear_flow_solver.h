////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
#define TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H

#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>
#include <gmm.h>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMapNetwork;

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

/*! @brief Saves the material constants on a vessel. */
struct VesselParameters {
  VesselParameters(double G0, double A0, double rho)
      : G0(G0), A0(A0), rho(rho) {}

  double G0;
  double A0;
  double rho;
};

inline double calculate_c0(double G0, double rho, double A0) {
  return std::sqrt(G0 / (2 * rho));
}

inline double calculate_Q(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return 0.5 * (w2 - w1) * std::pow((w1 + w2) / (8 * c0), 4) * A0;
}

inline double calculate_diff_Q_w1(double w1, double w2, double G0, double rho, double A0)
{
  const double c0 = calculate_c0(G0, rho, A0);
  return - 0.5 * std::pow((w1+w2)/(8*c0), 4) * A0 + 2 * A0 * (w2-w1)*std::pow(w2+w1, 3)/std::pow(8*c0, 4);
}

inline double calculate_diff_Q_w2(double w1, double w2, double G0, double rho, double A0)
{
  const double c0 = calculate_c0(G0, rho, A0);
  return + 0.5 * std::pow((w1+w2)/(8*c0), 4) * A0 + 2 * A0 * (w2-w1)*std::pow(w2+w1, 3)/std::pow(8*c0, 4);
}

inline double calculate_p(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return 0.5 * rho * std::pow(0.5*(w2 - w1), 2) + G0 * (std::pow((w1+w2)/(8*c0), 2) - 1) ;
}

inline double calculate_diff_p_w1(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return - 0.25 * rho * (w2-w1)  + 2 * G0 * (w1 + w2) / std::pow((8*c0), 2);
}

inline double calculate_diff_p_w2(double w1, double w2, double G0, double rho, double A0) {
  const double c0 = calculate_c0(G0, rho, A0);
  return + 0.25 * rho * (w2-w1)  + 2 * G0 * (w1 + w2) / std::pow((8*c0), 2);
}

inline void bifurcation_equation(
  double w1_a, double w2_a, const VesselParameters &p_a, bool in_a, // flow vessel a
  double w1_b, double w2_b, const VesselParameters &p_b, bool in_b, // flow vessel b
  double w1_c, double w2_c, const VesselParameters &p_c, bool in_c, // flow vessel c
  std::vector<double> &out) {

  assert(out.size() == 3);
  double eq_Q = 0;
  eq_Q += (in_a ? -1 : +1) * calculate_Q(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);
  eq_Q += (in_b ? -1 : +1) * calculate_Q(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);
  eq_Q += (in_c ? -1 : +1) * calculate_Q(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);
  out[0] = eq_Q;
  double eq_p1 = calculate_p(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0) - calculate_p(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);
  double eq_p2 = calculate_p(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0) - calculate_p(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);
  out[0] = eq_Q;
  out[1] = eq_p1;
  out[2] = eq_p2;
}

inline void bifurcation_equation_jacobian(
  double w1_a, double w2_a, const VesselParameters &p_a, bool in_a, // flow vessel a
  double w1_b, double w2_b, const VesselParameters &p_b, bool in_b, // flow vessel b
  double w1_c, double w2_c, const VesselParameters &p_c, bool in_c, // flow vessel c
  gmm::row_matrix<gmm::wsvector<double>>& J) {

  assert(J.ncols() == 3 && J.nrows() == 3);

  // line 0
  if (in_a)
    J[0][0] = - calculate_diff_Q_w1(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);
  else
    J[0][0] = + calculate_diff_Q_w2(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);

  if (in_b)
    J[0][1] = - calculate_diff_Q_w1(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);
  else
    J[0][1] = + calculate_diff_Q_w2(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);

  if (in_c)
    J[0][2] = - calculate_diff_Q_w1(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);
  else
    J[0][2] = + calculate_diff_Q_w2(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);

  // line 1
  if (in_a)
    J[1][0] = + calculate_diff_p_w1(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);
  else
    J[1][0] = + calculate_diff_p_w2(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);

  if (in_b)
    J[1][1] = - calculate_diff_p_w1(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);
  else
    J[1][1] = - calculate_diff_p_w2(w1_b, w2_b, p_b.G0, p_b.rho, p_b.A0);

  // line 2
  if (in_a)
    J[2][0] = + calculate_diff_p_w1(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);
  else
    J[2][0] = + calculate_diff_p_w2(w1_a, w2_a, p_a.G0, p_a.rho, p_a.A0);

  if (in_c)
    J[2][2] = - calculate_diff_p_w1(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);
  else
    J[2][2] = - calculate_diff_p_w2(w1_c, w2_c, p_c.G0, p_c.rho, p_c.A0);
}

inline void extract_vector(double& w1_a, double& w2_a, bool in_a, double& w1_b, double& w2_b, bool in_b, double& w1_c, double& w2_c, bool in_c, const std::vector< double > & vec) {
  if (in_a)
    w1_a = vec[0];
  else
    w2_a = vec[0];

  if (in_b)
    w1_b = vec[1];
  else
    w2_b = vec[1];

  if (in_c)
    w1_c = vec[2];
  else
    w2_c = vec[2];
}

inline void fill_vector(double w1_a, double w2_a, bool in_a, double w1_b, double w2_b, bool in_b, double w1_c, double w2_c, bool in_c, std::vector< double > & vec) {
  if (in_a)
    vec[0] = w1_a;
  else
    vec[0] = w2_a;

  if (in_b)
    vec[1] = w1_b;
  else
    vec[1] = w2_b;

  if (in_c)
    vec[2] = w1_c;
  else
    vec[2] = w2_c;
}

/*! @brief Evaluates the back propagating wave. */
inline double calculate_W1_value(const double Q, const double A, const double G0, const double rho, const double A0) {
  return -Q / A + 4 * std::sqrt(G0 / (2 * rho)) * std::pow(A / A0, 1. / 4.);
}

/*! @brief Evaluates the forward propagating wave. */
inline double calculate_W2_value(double Q, double A, double G0, double rho, double A0) {
  return +Q / A + 4 * std::sqrt(G0 / (2 * rho)) * std::pow(A / A0, 1. / 4.);
}

inline double convert_w1w2_to_QA(double w1, double w2, const VesselParameters & p, double& Q, double& A)
{
  const double c0 = calculate_c0(p.G0, p.rho, p.A0);
  A = std::pow((w1+w2)/(8*c0), 4) * p.A0;
  Q = 0.5*(w2 - w1)*A;
}

inline std::size_t solve_at_bifurcation(
  double Q_a, double A_a, const VesselParameters &p_a, bool in_a, // flow vessel a
  double Q_b, double A_b, const VesselParameters &p_b, bool in_b, // flow vessel b
  double Q_c, double A_c, const VesselParameters &p_c, bool in_c, // flow vessel c
  double& Q_a_up, double& A_a_up,
  double& Q_b_up, double& A_b_up,
  double& Q_c_up, double& A_c_up)
{
  const std::size_t max_iter = 1000;

  // calculate w1, w2 at all the bifurcations
  double w1_a = calculate_W1_value(Q_a, A_a, p_a.G0, p_a.rho, p_a.A0);
  double w2_a = calculate_W2_value(Q_a, A_a, p_a.G0, p_a.rho, p_a.A0);
  double w1_b = calculate_W1_value(Q_b, A_b, p_b.G0, p_b.rho, p_b.A0);
  double w2_b = calculate_W2_value(Q_b, A_b, p_b.G0, p_b.rho, p_b.A0);
  double w1_c = calculate_W1_value(Q_c, A_c, p_c.G0, p_c.rho, p_c.A0);
  double w2_c = calculate_W2_value(Q_c, A_c, p_c.G0, p_c.rho, p_c.A0);

  // solve for w1, w2 with a newton iteration
  gmm::row_matrix<gmm::wsvector<double>> Jac(3, 3);
  std::vector< double > x(3, 0);
  std::vector< double > f(3, 0);
  std::vector< double > delta_x(3, 0);

  fill_vector(w1_a, w2_a, in_a, w1_b, w2_b, in_b, w1_c, w2_c, in_c, x);
  std::size_t it=0;
  for (; it<max_iter; it+=1)
  {
    bifurcation_equation_jacobian(w1_a, w2_a, p_a, in_a, w1_b, w2_b, p_b, in_b, w1_c, w2_c, p_c, in_c, Jac);
    bifurcation_equation(w1_a, w2_a, p_a, in_a, w1_b, w2_b, p_b, in_b, w1_c, w2_c, p_c, in_c, f);
    gmm::lu_solve(Jac, delta_x, f);
    gmm::add(x, gmm::scaled(delta_x, -1), x);
    extract_vector(w1_a, w2_a, in_a, w1_b, w2_b, in_b, w1_c, w2_c, in_c, x);

    if (gmm::vect_norm2(delta_x) < 1e-14 || gmm::vect_norm2(delta_x) < 1e-14 * gmm::vect_norm2(x))
      break;
  }

  // copy back into the input variables
  convert_w1w2_to_QA(w1_a, w2_a, p_a, Q_a_up, A_a_up);
  convert_w1w2_to_QA(w1_b, w2_b, p_b, Q_b_up, A_b_up);
  convert_w1w2_to_QA(w1_c, w2_c, p_c, Q_c_up, A_c_up);

  return it;
}

/*! @brief Solves the system of equation
 *           W1(Q_up, A_up) = W1
 *           W2(Q_up, A_up) = W2
 *         for (Q_up, Q_up).
 */
inline void solve_W12(double &Q_up, double &A_up, const double W1, const double W2, const double G0, const double rho, const double A0) {
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
inline double heart_beat_inflow(double t) {
  const double t_period = 1.0;
  const double t_systole = 0.3;
  const double t_in_period = t - std::floor(t / t_period);
  if (t_in_period < t_systole) {
    return 485. * std::sin(M_PI * t_in_period / t_systole);
  } else if (t_in_period <= t_period + 1e-14) {
    return 0.;
  } else {
    throw std::runtime_error("unreachable code");
  }
}

inline double solve_for_W2(const double W1, const double W2_init, const double Q_star, const double A_0, const double c_0) {
  double W2 = W2_init;
  double W2_prev = W2_init;
  const auto f = [=](double W2) { return (W2 - W1) / 2. * std::pow((W1 + W2) / (8 * c_0), 4) * A_0 - Q_star; };
  const auto f_prime = [=](double W2) { return (0.000610352 * A_0 * (W2 - 0.6 * W1) * pow((W2 + W1), 3)) / pow(c_0, 4); };
  for (std::size_t it = 0; it < 100; it += 1) {
    W2 += (-1) * f(W2) / f_prime(W2);
    if (std::abs(W2 - W2_prev) < 1e-16 || std::abs(W2 - W2_prev) < 1e-8 * W2_prev)
      break;
    W2_prev = W2;
  }
  return W2;
}

inline double assemble_in_flow(const double Q_l, const double A_l, const double Q_star, const double G0, const double rho, const double A0) {
  const double c0 = sqrt(G0 / (2 * rho));
  const double W1 = calculate_W1_value(Q_l, A_l, G0, rho, A0);
  const double W2_init = calculate_W2_value(Q_l, A_l, G0, rho, A0);
  const double W2 = solve_for_W2(W1, W2_init, Q_star, A0, c0);
  const double A = A0 * std::pow(1. / (8 * c0) * (W2 + W1), 4);
  return A;
}

/*! @brief Assembles the inverse mass. WARNING: Assumes legendre basis! */
void assemble_inverse_mass(const GraphStorage &graph, const DofMapNetwork &dof_map, std::vector<double> &inv_mass);

/*! @brief Interpolates a constant value. WARNING: Assumes legendre basis! */
void interpolate_constant(const GraphStorage &graph, const DofMapNetwork &dof_map, double value, std::size_t component, std::vector<double> &result);

class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph);

  void solve();

  double get_time() const;

  void set_tau(double tau);

  double get_solution_on_vertices(std::vector<double> &Q_values, std::vector<double> &A_values) const;

private:
  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMapNetwork> d_dof_map;

  /*! @brief Current time step size. */
  double d_tau;

  /*! @brief The current time. */
  double d_t_now;

  std::function<double(double)> d_inflow_value_function;

  /*! @brief The solution at the current time step. */
  std::vector<double> d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector<double> d_u_prev;

  /*! @brief The upwinded Q-flow. d_Q_up[vertex_id] contains the respective flux. */
  std::vector<double> d_Q_up;

  /*! @brief The upwinded A-flow. d_A_up[vertex_id] contains the respective flux. */
  std::vector<double> d_A_up;

  /*! @brief Our current rhs vector before multiplying it with the inverse mass. */
  std::vector<double> d_rhs;

  /*! @brief Our current inverse mass vector, defining the diagonal inverse mass matrix. */
  std::vector<double> d_inverse_mass;

  // All the parameters, which are now kept constant, but should change for vessel types.
  double d_G0;
  double d_rho;
  double d_A0;
  double d_mu;
  double d_gamma;

  /*! @brief Recalculates the current fluxes from the previous time step. */
  void calculate_fluxes();

  /*! @brief Assembles from the fluxes and the previous values a new right hand side function. */
  void calculate_rhs();

  /*! @brief Applies the inverse mass to our right hand side and saves the result to u_now. */
  void apply_inverse_mass();
};

} // namespace macrocirculation

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
