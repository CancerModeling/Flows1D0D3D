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

  /*! @brief The upwinded Q-flow from the left boundary of a macroedge,
   *         ordered such that d_Q_up_el[edge_id] contains the respective flux.
   * */
  std::vector<double> d_Q_up_el;

  /*! @brief The upwinded Q-flow from the right boundary of a macroedge,
   *         ordered such that d_Q_up_er[edge_id] contains the respective flux.
   * */
  std::vector<double> d_Q_up_er;

  /*! @brief The upwinded A-flow from the left boundary of a macroedge,
   *         ordered such that d_A_up_el[edge_id] contains the respective flux.
   * */
  std::vector<double> d_A_up_el;

  /*! @brief The upwinded A-flow from the right boundary of a macroedge,
   *         ordered such that d_Q_up_er[edge_id] contains the respective flux.
   * */
  std::vector<double> d_A_up_er;

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
