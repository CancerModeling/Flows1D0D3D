////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.h"
#include <cassert>

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"
#include "interpolate_to_vertices.hpp"
#include "right_hand_side_evaluator.hpp"
#include "time_integrators.hpp"
#include "vessel_data_storage.hpp"

namespace macrocirculation {

namespace lm = libMesh;

void interpolate_constant(const GraphStorage &graph, const DofMapNetwork &dof_map, double value, std::size_t component, std::vector<double> &result) {
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    dof_map.dof_indices(*edge, dof_indices, component);
    result[dof_indices[0]] = value;
  }
}

void set_to_A0(const GraphStorage &graph, const DofMapNetwork &dof_map, const VesselDataStorage& vessel_data, std::vector<double> &result)
{
  std::vector<std::size_t> dof_indices;
  for (const auto &e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    const auto& data = vessel_data.get_parameters(*edge);
    // set Q
    dof_map.dof_indices(*edge, dof_indices, 0);
    result[dof_indices[0]] = 0;
    // set A
    dof_map.dof_indices(*edge, dof_indices, 1);
    result[dof_indices[0]] = data.A0;
  }
}

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data)
    : d_graph(std::move(graph)),
      d_vessel_data(std::move(vessel_data)),
      d_dof_map(std::make_shared<SimpleDofMapNetwork>(2, degree + 1, d_graph->num_edges())),
      d_right_hand_side_evaluator(std::make_shared<RightHandSideEvaluator<degree>>(d_graph, d_vessel_data, d_dof_map)),
      d_time_integrator(std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof())),
      d_tau(2.5e-4 / 4),
      d_t_now(0),
      d_u_now(d_dof_map->num_dof()),
      d_u_prev(d_dof_map->num_dof())
{
  // set A constant to A0
  set_to_A0(*d_graph, *d_dof_map, *d_vessel_data, d_u_prev);
  set_to_A0(*d_graph, *d_dof_map, *d_vessel_data, d_u_now);
}

// we need the destructor here, to use unique_ptrs with forward declared classes.
ExplicitNonlinearFlowSolver::~ExplicitNonlinearFlowSolver() = default;

void ExplicitNonlinearFlowSolver::solve() {
  d_u_prev = d_u_now;
  d_t_now += d_tau;
  const double t_prev = d_t_now - d_tau;
  d_time_integrator->apply<degree>(d_u_prev, t_prev, d_tau, *d_right_hand_side_evaluator, d_u_now);
}

double ExplicitNonlinearFlowSolver::get_time() const { return d_t_now; }

void ExplicitNonlinearFlowSolver::set_tau(double tau) { d_tau = tau; }

void ExplicitNonlinearFlowSolver::use_explicit_euler_method()
{
  d_time_integrator = std::make_unique<TimeIntegrator>(create_explicit_euler(), d_dof_map->num_dof());
}

void ExplicitNonlinearFlowSolver::use_ssp_method()
{
  d_time_integrator = std::make_unique<TimeIntegrator>(create_ssp_method(), d_dof_map->num_dof());
}

void ExplicitNonlinearFlowSolver::get_solution_on_vertices(std::vector<double> &Q_values, std::vector<double> &A_values) const {
  assert(Q_values.size() == d_graph->num_edges() * 2);
  assert(A_values.size() == d_graph->num_edges() * 2);

  FETypeNetwork<degree> fe(create_trapezoidal_rule());

  interpolate_to_vertices(*d_graph, *d_dof_map, fe, 0, d_u_now, Q_values);
  interpolate_to_vertices(*d_graph, *d_dof_map, fe, 1, d_u_now, A_values);
}

} // namespace macrocirculation