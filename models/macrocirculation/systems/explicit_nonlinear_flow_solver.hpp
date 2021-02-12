////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_HPP
#define TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_HPP

#include <cassert>
#include <cmath>
#include <functional>
#include <gmm.h>
#include <memory>
#include <vector>

namespace macrocirculation {

constexpr std::size_t degree = 1;

// forward declarations
class GraphStorage;
class VesselDataStorage;
class DofMapNetwork;
template<std::size_t degree>
class RightHandSideEvaluator;
class TimeIntegrator;

/*! @brief Interpolates a constant value. WARNING: Assumes legendre basis! */
void interpolate_constant(const GraphStorage &graph, const DofMapNetwork &dof_map, double value, std::size_t component, std::vector<double> &result);

/*! @brief Sets the given function to A=A0 and Q=0. WARNING: Assumes legendre basis! */
void set_to_A0(const GraphStorage &graph, const DofMapNetwork &dof_map, const VesselDataStorage&, std::vector<double> &result);

class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data_storage);
  ~ExplicitNonlinearFlowSolver();

  void solve();

  double get_time() const;

  void set_tau(double tau);

  void get_solution_on_vertices(std::vector<double> &Q_values, std::vector<double> &A_values) const;

  void get_total_pressure_on_vertices(std::vector<double> &p_values) const;

  void get_static_pressure_on_vertices(std::vector<double> &p_values) const;

  /*! @brief Configures the explicit euler method as the time integrator. */
  void use_explicit_euler_method();

  /*! @brief Configures a 3rd order RKM as the time integrator. */
  void use_ssp_method();

private:
  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief Stores the physical data for a given vessel. */
  std::shared_ptr<VesselDataStorage> d_vessel_data;

  /*! @brief The dof map for our domain. */
  std::shared_ptr<DofMapNetwork> d_dof_map;

  /*! @brief Utility class for evaluating the right hand side, to allow different explicit schemes. */
  std::shared_ptr<RightHandSideEvaluator<degree>> d_right_hand_side_evaluator;

  /*! @brief Explicit time integrator to move the solution forwards in time. */
  std::unique_ptr<TimeIntegrator> d_time_integrator;

  /*! @brief Current time step size. */
  double d_tau;

  /*! @brief The current time. */
  double d_t_now;

  /*! @brief The solution at the current time step. */
  std::vector<double> d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector<double> d_u_prev;

  /*! @brief The solution at the previous time step. */
  std::vector< std::vector<double> > d_k;
};

} // namespace macrocirculation

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_HPP
