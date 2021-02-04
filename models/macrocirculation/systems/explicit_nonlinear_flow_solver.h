////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
#define TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H

#include "libmesh/libmesh.h"

#include <vector>
#include <cmath>
#include <memory>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMapNetwork;

/*! @brief Evaluates the back propagating wave. */
inline const double get_W1_value(const double Q, const double A, const double G0, const double rho, const double A0 )
{
  return - Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
}

/*! @brief Evaluates the forward propagating wave. */
inline double get_W2_value(double Q, double A, double G0, double rho, double A0 )
{
  return + Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
}

/*! @brief Calculates the Jacobian of the (W1, W2) vector valued function. */
template < typename Mat >
inline const void get_W12_jacobian(const double Q, const double A, const double G0, const double rho, const double A0, Mat& mat)
{
  mat[0, 0] = - 1./A;
  mat[1, 0] = + 1./A;
  const double a = Q/std::pow(A, 2);
  const double b = std::sqrt(G0/rho) / ( std::sqrt(2) * A0 * std::pow(A/A0, 3./4.));
  mat[0, 1] =  a + b;
  mat[1, 1] = -a + b;
}

/*! @brief Solves the system of equation
 *           W1(Q_up, A_up) = W1
 *           W2(Q_up, A_up) = W2
 *         for (Q_up, Q_up) by applying a Newton iteration.
 */
template < typename Mat >
inline const void solve_W12(double & Q_up, double & A_up, const double W1, const double W2,  const double G0, const double rho, const double A0)
{
}

class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver( std::shared_ptr<GraphStorage> graph );

  void solve();

private:
  /*! @brief Current time step size. */
  double d_tau;

  /*! @brief The current time. */
  double d_t_now;

  /*! @brief The solution at the current time step. */
  std::vector< double > d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector< double > d_u_prev;

  /*! @brief The upwinded flows. */
  std::vector< double > d_u_up;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr< GraphStorage > d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr< DofMapNetwork > d_dof_map;

  /*! @brief Recalculates the current fluxes from the previous time step. */
  void calculate_fluxes();

  // All the parameters, which are now kept constant, but should change for vessel types.
  double d_G0;
  double d_rho;
  double d_A0;
};

}

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
