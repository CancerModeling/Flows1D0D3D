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

/*! @brief Evaluates the back propagating wave. */
inline double get_W1_value(double Q, double A, double G0, double rho, double A0 )
{
  return - Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
}

/*! @brief Evaluates the forward propagating wave. */
inline double get_W2_value(double Q, double A, double G0, double rho, double A0 )
{
  return + Q/A + 4 * std::sqrt(G0/(2*rho)) * std::pow(A/A0, 1./4.);
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

  /*! @brief Recalculates the current fluxes from the previous time step. */
  void calculate_fluxes();

  // All the parameters, which are now kept constant, but should change for vessel types.
  double d_G0;
  double d_rho;
  double d_A0;
};

}

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_H
