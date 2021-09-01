////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_NONLINEAR_FLOW_UPWIND_EVALUATOR_HPP
#define TUMORMODELS_NONLINEAR_FLOW_UPWIND_EVALUATOR_HPP

#include <mpi.h>
#include <memory>
#include <vector>

#include "edge_boundary_evaluator.hpp"

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class Vertex;
class Edge;

/*! @brief Class for calculating the currently upwinded values for the flow (Q, A).
 *
 *  This is in a separate class from the right-hand-side evaluator, since the transport equation also needs these flows.
 *  And since transport and flow can be decoupled it makes sense to keep them apart.
 */
class NonlinearFlowUpwindEvaluator {
public:
  NonlinearFlowUpwindEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map);

  void init(double t, const std::vector<double> &u_prev);

  /*! @brief Recalculates the current fluxes from the given time.
   *
   * @param t       The current time for the inflow boundary conditions.
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void get_fluxes_on_macro_edge(double t, const Edge &edge, const std::vector<double> &u_prev, std::vector<double> &Q_up, std::vector<double> &A_up) const;

  /*! @brief Recalculates the current flux values of (Q, A) at the given nfurcation.
   *
   * @param t       The current time for the inflow boundary conditions.
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &Q_up, std::vector<double> &A_up) const;

private:
  /*! @brief Calculates the fluxes at nfurcations for the given time step at the macro edge boundaries.
   *
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_nfurcation_fluxes(const std::vector<double> &u_prev);

  /*! @brief Calculates the in and outflow fluxes for the given timestep at the macro edge boundaries.
   *
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_inout_fluxes(double t, const std::vector<double> &u_prev);

private:
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMap> d_dof_map;

  /*! @brief Stores the values of the A-component at the macro-edge boundaries. */
  EdgeBoundaryEvaluator d_A_boundary_evaluator;

  /*! @brief Stores the values of the Q-component at the macro-edge boundaries. */
  EdgeBoundaryEvaluator d_Q_boundary_evaluator;

  /*! @brief Contains the flux-values of Q at the left boundary point of the given macro-edge.
    *         The ith vector entry corresponds to the ith macro-edge id.
    */
  std::vector<double> d_Q_macro_edge_flux_l;
  std::vector<double> d_Q_macro_edge_flux_r;
  std::vector<double> d_A_macro_edge_flux_l;
  std::vector<double> d_A_macro_edge_flux_r;

  double d_current_t;
};

} // namespace macrocirculation

#endif