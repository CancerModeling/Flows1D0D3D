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
#include <mpi.h>
#include <vector>

namespace macrocirculation {

constexpr std::size_t degree = 2;

// forward declarations
class GraphStorage;
class DofMap;
class RightHandSideEvaluator;
class TimeIntegrator;
class Vertex;

/*! @brief Interpolates a constant value. WARNING: Assumes legendre basis! */
void interpolate_constant(MPI_Comm comm,
                          const GraphStorage &graph,
                          const DofMap &dof_map,
                          double value,
                          std::size_t component,
                          std::vector<double> &result);

/*! @brief Sets the given function to A=A0 and Q=0. WARNING: Assumes legendre basis! */
void set_to_A0(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, std::vector<double> &result);

template<std::size_t degree>
class ExplicitNonlinearFlowSolver {
public:
  explicit ExplicitNonlinearFlowSolver(MPI_Comm comm,
                                       std::shared_ptr<GraphStorage> graph,
                                       std::shared_ptr<DofMap> dof_map);
  ~ExplicitNonlinearFlowSolver();

  void solve();

  double get_time() const;

  void set_tau(double tau);

  /*! @brief Configures the explicit euler method as the time integrator. */
  void use_explicit_euler_method();

  /*! @brief Configures a 3rd order RKM as the time integrator. */
  void use_ssp_method();

  RightHandSideEvaluator &get_rhs_evaluator();

  DofMap &get_dof_map();

  std::vector<double> &get_solution();

  // TODO: Move this somewhere else
  /*! @brief Calculates the flow Q pointing towards the vertex v. */
  [[nodiscard]] double get_flow_at_vessel_tip(const Vertex& v) const;

private:
  /*! @brief The mpi communicator. */
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain. */
  std::shared_ptr<DofMap> d_dof_map;

  /*! @brief Utility class for evaluating the right hand side, to allow different explicit schemes. */
  std::shared_ptr<RightHandSideEvaluator> d_right_hand_side_evaluator;

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
  std::vector<std::vector<double>> d_k;
};

} // namespace macrocirculation

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_HPP
