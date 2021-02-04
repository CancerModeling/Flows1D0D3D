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

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMapNetwork;

/*! @brief Calculates Eq. (2.8) from "Multi-scale modeling of flow and transport processes in arterial networks and tissue".
 *
 * @param h0 The thickness of the blood vessel.
 * @param E Young's modulus.
 * @param nu Poisson ratio.
 * @param A0 The vessel area
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

class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver( std::shared_ptr<GraphStorage> graph );

  void solve();

private:
  /*! @brief The current domain for solving the equation. */
  std::shared_ptr< GraphStorage > d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr< DofMapNetwork > d_dof_map;

  /*! @brief Current time step size. */
  double d_tau;

  /*! @brief The current time. */
  double d_t_now;

  /*! @brief The solution at the current time step. */
  std::vector< double > d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector< double > d_u_prev;

  /*! @brief The upwinded Q-flow. d_Q_up[vertex_id] contains the respective flux. */
  std::vector< double > d_Q_up;

  /*! @brief The upwinded A-flow. d_A_up[vertex_id] contains the respective flux. */
  std::vector< double > d_A_up;

  /*! @brief Recalculates the current fluxes from the previous time step. */
  void calculate_fluxes();

  // All the parameters, which are now kept constant, but should change for vessel types.
  double d_G0;
  double d_rho;
  double d_A0;
};

}

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
