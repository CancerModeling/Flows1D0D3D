////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TRANSPORT_UPWIND_EVALUATOR_HPP
#define TUMORMODELS_TRANSPORT_UPWIND_EVALUATOR_HPP

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
class UpwindProvider;
class PetscVec;

/*! @brief Class for calculating the currently upwinded values for the transport. */
class TransportUpwindEvaluator {
public:
  TransportUpwindEvaluator(MPI_Comm comm,
                           std::shared_ptr<GraphStorage> graph,
                           std::shared_ptr<DofMap> dof_map,
                           std::shared_ptr<UpwindProvider> flow_upwind_provider);

  void init(double t, const PetscVec &u_prev);

  void get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &gamma_up) const;

  bool upwind_is_implemented(const Vertex& v) const;

  void set_inflow_function(std::function< double(double) > inflow_function);

private:
  /*! @brief Calculates the fluxes at nfurcations for the given time step at the macro edge boundaries.
   *
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_nfurcation_fluxes(double t, const PetscVec &u_prev);

  /*! @brief Calculates the in and outflow fluxes for the given time step at the macro edge boundaries.
   *
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_inout_fluxes(double t, const PetscVec &u_prev);

private:
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMap> d_dof_map;

  /*! @brief Stores the values at the macro-edge boundaries. */
  EdgeBoundaryEvaluator d_gamma_boundary_evaluator;

  std::shared_ptr< UpwindProvider > d_flow_upwind_provider;

  std::function< double(double) > d_inflow_function = [](double t){ return 1.; };

  /*! @brief Contains the flux-values of Q at the left boundary point of the given macro-edge.
    *         The ith vector entry corresponds to the ith macro-edge id.
    */
  std::vector<double> d_gamma_macro_edge_flux_l;
  std::vector<double> d_gamma_macro_edge_flux_r;

  double d_current_t;
};

} // namespace macrocirculation

#endif