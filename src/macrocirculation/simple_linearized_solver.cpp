////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2022 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "dof_map.hpp"
#include "graph_partitioner.hpp"
#include "graph_pvd_writer.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "simple_linearized_solver.hpp"
#include "vessel_formulas.hpp"
#include "embedded_graph_reader.hpp"
#include "petsc.h"

namespace macrocirculation {

SimpleLinearizedSolver::SimpleLinearizedSolver(const std::string& filepath)
// TODO: Manage parallelism better
  : comm(PETSC_COMM_SELF), tau(1e-3), t(0), writer(std::make_shared<GraphPVDWriter>(comm, "./", "simple_linearized_solver")) {

  graph = std::make_shared<GraphStorage>();

  EmbeddedGraphReader graph_reader;
  graph_reader.append(filepath, *graph);

  set_tau(tau);
}

void SimpleLinearizedSolver::set_tau(double ptau) {
  tau = ptau;

  const size_t degree = 2;

  v_coupling_1_outer = graph->find_vertex_by_name("coupling_1_outer");
  v_coupling_1_inner = graph->find_vertex_by_name("coupling_1_inner");
  v_coupling_2_inner = graph->find_vertex_by_name("coupling_2_inner");
  v_coupling_2_outer = graph->find_vertex_by_name("coupling_2_outer");

  edge0 = graph->find_edge_by_name("vessel_coupling_1");
  edge1 = graph->find_edge_by_name("vessel_coupling_2");

  if (graph->has_named_vertex("inflow"))
  {
    auto v0 = graph->find_vertex_by_name("inflow");
    v0->set_to_inflow_with_fixed_pressure([](double t) { return 2 * std::abs(std::sin(M_PI * t)); });
  }

  if (graph->has_named_vertex("outflow"))
  {
    auto v3 = graph->find_vertex_by_name("outflow");
    v3->set_to_free_outflow();
  }

  graph->finalize_bcs();

  naive_mesh_partitioner(*graph, comm);

  dof_map = std::make_shared<DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(comm, *graph, 2, degree, true);

  solver = std::make_shared<ImplicitLinearFlowSolver>(comm, graph, dof_map, degree);
  solver->setup(tau);
}

void SimpleLinearizedSolver::solve() {
  t += tau;
  solver->solve(tau, t);
}

SimpleLinearizedSolver::Result SimpleLinearizedSolver::get_result(const Vertex &vertex, const Edge &edge) {
  double p, q;
  solver->get_1d_pq_values_at_vertex(vertex, edge, p, q);
  p *= 1e3; // mixed units to cgs
  const auto &phys_data = edge.get_physical_data();
  const double a = nonlinear::get_A_from_p(p, phys_data.G0, phys_data.A0);
  return {a, p, q};
}

SimpleLinearizedSolver::Result SimpleLinearizedSolver::get_result(Outlet outlet) {
  if (outlet == Outlet::in) {
    return get_result(*v_coupling_1_inner, *edge0);
  } else if (outlet == Outlet::out) {
    return get_result(*v_coupling_2_inner, *edge1);
  } else {
    throw std::runtime_error("unknown vertex value");
  }
}

void SimpleLinearizedSolver::write() {
  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;
  interpolate_to_vertices(comm, *graph, *dof_map, solver->p_component, solver->get_solution(), points, p_vertex_values);
  interpolate_to_vertices(comm, *graph, *dof_map, solver->q_component, solver->get_solution(), points, q_vertex_values);

  writer->set_points(points);
  writer->add_vertex_data("p", p_vertex_values);
  writer->add_vertex_data("q", q_vertex_values);
  writer->write(t);
}

} // namespace macrocirculation