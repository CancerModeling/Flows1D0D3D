////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TRANSPORT_H
#define TUMORMODELS_TRANSPORT_H

#include <cmath>
#include <memory>
#include <mpi.h>
#include <utility>
#include <vector>

#include "flow_upwind_evaluator.hpp"

namespace macrocirculation {

// forward declaration
class GraphStorage;
class DofMap;
class Edge;

class Transport {
public:
  Transport(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map_flow, std::shared_ptr<DofMap> dof_map_transport);

  void evaluate_macro_edge_boundary_values(const std::vector<double> &u_prev, const std::vector<double> &gamma_prev);

  std::vector<double> &get_solution();

  void solve(double t, double dt, const std::vector<double> &u_prev);

  void calculate_fluxes_on_macro_edge(double t,
                                      const Edge &edge,
                                      const std::vector<double> &u_prev,
                                      const std::vector<double> &gamma_prev,
                                      std::vector<double> &gamma_fluxes_edge);

  void assemble_rhs(double t, const std::vector<double> &u_prev, const std::vector<double> &gamma_prev, std::vector<double> &rhs);

private:
  void calculate_fluxes_at_nfurcations(double t, const std::vector<double> &u_prev);

  void apply_inverse_mass();

private:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dof_map_flow;
  std::shared_ptr<DofMap> d_dof_map_transport;

  FlowUpwindEvaluator d_flow_upwind_evaluator;

  Communicator d_edge_boundary_communicator;

  std::vector<double> d_gamma_macro_edge_boundary_value;

  std::vector<double> d_gamma_flux_l;
  std::vector<double> d_gamma_flux_r;

  std::vector<double> d_rhs;
  std::vector<double> d_solution;
};

} // namespace macrocirculation

#endif //TUMORMODELS_TRANSPORT_H
