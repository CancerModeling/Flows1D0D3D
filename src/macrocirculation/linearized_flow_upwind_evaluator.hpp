////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_LINEARIZED_FLOW_UPWIND_EVALUATOR_HPP
#define TUMORMODELS_LINEARIZED_FLOW_UPWIND_EVALUATOR_HPP

#include "edge_boundary_evaluator.hpp"

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class DofMap;
class PetscVec;
class Edge;
class Vertex;

class LinearizedFlowUpwindEvaluator {
public:
  LinearizedFlowUpwindEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map);

  void init(double t, const std::vector<double> &u_prev);

  void init(double t, const PetscVec &u_prev);

  void get_fluxes_on_macro_edge(double t, const Edge &edge, const std::vector<double> &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const;

  void get_fluxes_on_macro_edge(double t, const Edge &edge, const PetscVec &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const;

  void get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &p_up, std::vector<double> &q_up) const;

private:
  template<typename VectorType>
  void init_generic(double t, const VectorType &u_prev);

  template<typename VectorType>
  void get_fluxes_on_macro_edge_generic(double t, const Edge &edge, const VectorType &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const;

  template<typename VectorType>
  void calculate_nfurcation_fluxes(const VectorType &u_prev);

  template<typename VectorType>
  void calculate_inout_fluxes(double t, const VectorType &u_prev);

private:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<DofMap> d_dof_map;

  EdgeBoundaryEvaluator d_p_boundary_evaluator;

  EdgeBoundaryEvaluator d_q_boundary_evaluator;

  std::vector<double> d_p_macro_edge_flux_l;
  std::vector<double> d_p_macro_edge_flux_r;
  std::vector<double> d_q_macro_edge_flux_l;
  std::vector<double> d_q_macro_edge_flux_r;

  double d_current_t;
};

} // namespace macrocirculation

#endif //TUMORMODELS_LINEARIZED_FLOW_UPWIND_EVALUATOR_HPP
