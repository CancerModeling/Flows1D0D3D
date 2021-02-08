////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
#define TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H

#include <vector>
#include <cmath>
#include <memory>
#include <functional>

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
inline double calculate_G0(double h0, double E, double nu, double A0)
{
  return std::sqrt(M_PI)*h0*E / ((1-nu*nu)*std::sqrt(A0));
}

/*! @brief Evaluates the back propagating wave. */
inline double calculate_W1_value(const double Q, const double A, const double G0, const double rho, const double A0 )
{
  return - Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
}

/*! @brief Evaluates the forward propagating wave. */
inline double calculate_W2_value(double Q, double A, double G0, double rho, double A0 )
{
  return + Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
}

/*! @brief Solves the system of equation
 *           W1(Q_up, A_up) = W1
 *           W2(Q_up, A_up) = W2
 *         for (Q_up, Q_up).
 */
inline void solve_W12(double & Q_up, double & A_up, const double W1, const double W2,  const double G0, const double rho, const double A0)
{
  const double in = 1./8. * std::sqrt(2. * rho / G0) * (W2 + W1);
  A_up = A0 * std::pow(in, 4);
  Q_up = A_up/2. * (W2 - W1);
}

/*! @brief Inflow boundary condition modeling a heart beat.
 *         These values are from T. Koeppls doctoral thesis, Eq. (2.124)
 *
 * @param t Current time.
 * @return The current flow rate.
 */
inline double heart_beat_inflow(double t)
{
  const double t_period = 1.0;
  const double t_systole = 0.3;
  const double t_in_period = t - std::floor(t/t_period);
  if (t_in_period < t_systole)
  {
    return 485. * std::sin(M_PI * t_in_period/t_systole);
  }
  else if ( t_in_period <= t_period + 1e-14)
  {
    return 0.;
  }
  else
  {
    throw std::runtime_error("unreachable code");
  }
}

inline double solve_for_W2(const double W1, const double W2_init, const double Q_star, const double A_0, const double c_0)
{
  double W2 = W2_init;
  double W2_prev = W2_init;
  const auto f = [=](double W2){ return (W2 - W1) / 2. * std::pow((W1 + W2)/(8 * c_0), 4) * A_0 - Q_star; };
  const auto f_prime = [=](double W2){ return  (0.000610352*A_0*(W2 - 0.6*W1)*pow((W2 + W1), 3))/pow(c_0, 4); };
  for (std::size_t it=0; it<100; it+=1)
  {
    W2 += (-1) * f(W2)/f_prime(W2);
    if (std::abs(W2 - W2_prev) < 1e-16 || std::abs(W2 - W2_prev) < 1e-8 * W2_prev)
      break;
    W2_prev = W2;
  }
  return W2;
}

inline double assemble_in_flow(const double Q_l, const double A_l, const double Q_star, const double G0, const double rho, const double A0)
{
  const double c0 = sqrt(G0/(2*rho));
  const double W1 = calculate_W1_value(Q_l, A_l, G0, rho, A0);
  const double W2_init = calculate_W2_value(Q_l, A_l, G0, rho, A0);
  const double W2 = solve_for_W2(W1, W2_init, Q_star, A0, c0);
  const double A = A0*std::pow(1./(8*c0)*(W2+W1), 4);
  return A;
}

/*! @brief Assembles the inverse mass. WARNING: Assumes legendre basis! */
void assemble_inverse_mass(const GraphStorage & graph, const DofMapNetwork & dof_map, std::vector< double > & inv_mass);

/*! @brief Interpolates a constant value. WARNING: Assumes legendre basis! */
void interpolate_constant(const GraphStorage & graph, const DofMapNetwork & dof_map, double value, std::size_t component, std::vector< double > & result);

class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver( std::shared_ptr<GraphStorage> graph );

  void solve();

  double get_time() const;

  void set_tau(double tau);

  double get_solution_on_vertices(std::vector< double > & Q_values, std::vector< double > & A_values) const;

private:
  /*! @brief The current domain for solving the equation. */
  std::shared_ptr< GraphStorage > d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr< DofMapNetwork > d_dof_map;

  /*! @brief Current time step size. */
  double d_tau;

  /*! @brief The current time. */
  double d_t_now;

  std::function< double(double) > d_inflow_value_function;

  /*! @brief The solution at the current time step. */
  std::vector< double > d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector< double > d_u_prev;

  /*! @brief The upwinded Q-flow. d_Q_up[vertex_id] contains the respective flux. */
  std::vector< double > d_Q_up;

  /*! @brief The upwinded A-flow. d_A_up[vertex_id] contains the respective flux. */
  std::vector< double > d_A_up;

  /*! @brief Our current rhs vector before multiplying it with the inverse mass. */
  std::vector< double > d_rhs;

  /*! @brief Our current inverse mass vector, defining the diagonal inverse mass matrix. */
  std::vector< double > d_inverse_mass;

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

}

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
