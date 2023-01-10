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
#include "gmm_legacy_facade.hpp"
#include <memory>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class RightHandSideEvaluator;
class TimeIntegrator;
class Vertex;
class Edge;

struct Values0DModel {
  double p_c;
  double q;
};

/*! @brief Interpolates a constant value. WARNING: Assumes legendre basis! */
void interpolate_constant(MPI_Comm comm,
                          const GraphStorage &graph,
                          const DofMap &dof_map,
                          double value,
                          std::size_t component,
                          std::vector<double> &result);

/*! @brief Sets the given function to A=A0 and Q=0. WARNING: Assumes legendre basis! */
void set_to_A0(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, std::vector<double> &result);

/*! @brief An explicit solver for the nonlinear 1D flow equations. */
class ExplicitNonlinearFlowSolver {
public:
  /*! @brief Constructs an explicit nonlinear solver for the 1D flow equations.
   *
   * @param comm The parallel communicator for the solver.
   * @param graph The graph on which we solve the equations.
   * @param dof_map The dof map for the flow problem.
   * @param degree The degree of the finite element basis functions.
   */
  explicit ExplicitNonlinearFlowSolver(MPI_Comm comm,
                                       std::shared_ptr<GraphStorage> graph,
                                       std::shared_ptr<DofMap> dof_map,
                                       size_t degree);

  ~ExplicitNonlinearFlowSolver();

  /*! @brief The component index of the flow Q inside the dof map.
   *         This can be used inside a local dof map to get the dof indices for a single component. */
  static const size_t Q_component = 0;

  /*! @brief The component index of the area A inside the dof map.
   *         This can be used inside a local dof map to get the dof indices for a single component. */
  static const size_t A_component = 1;

  void solve(double tau, double t);

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

  void get_1d_AQ_values_at_vertex(const Vertex& v, double& A, double& Q) const;

  void get_1d_pq_values_at_vertex(const Vertex& v, double& p, double& q) const;

  [[nodiscard]] Values0DModel get_0D_values(const Vertex& v) const;

  /*! @brief Evaluates A and Q of the current solution on the edge e parametrized on [0, 1] at \f$ s \in [0,1] \f$. */
  void evaluate_1d_AQ_values(const Edge& e, double s, double& A, double& Q) const;

  /*! @brief Evaluates p and q of the current solution on the edge e parametrized on [0, 1] at \f$ s \in [0,1] \f$. */
  void evaluate_1d_pq_values(const Edge& e, double s, double& p, double& q) const;

  size_t get_degree() const;

private:
  /*! @brief The mpi communicator. */
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain. */
  std::shared_ptr<DofMap> d_dof_map;

  /*! @brief The degree of the finite element space. */
  size_t d_degree;

  /*! @brief Utility class for evaluating the right hand side, to allow different explicit schemes. */
  std::shared_ptr<RightHandSideEvaluator> d_right_hand_side_evaluator;

  /*! @brief Explicit time integrator to move the solution forwards in time. */
  std::unique_ptr<TimeIntegrator> d_time_integrator;

  /*! @brief The solution at the current time step. */
  std::vector<double> d_u_now;

  /*! @brief The solution at the previous time step. */
  std::vector<double> d_u_prev;

  /*! @brief The solution at the previous time step. */
  std::vector<std::vector<double>> d_k;
};

} // namespace macrocirculation

#endif //TUMORMODELS_EXPLICIT_NONLINEAR_FLOW_SOLVER_HPP
