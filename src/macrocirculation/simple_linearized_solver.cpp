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
#include "petsc.h"

namespace macrocirculation {

SimpleLinearizedSolver::SimpleLinearizedSolver()
    // TODO: Manage parallelism better
    : comm(PETSC_COMM_SELF), tau(1e-3), t(0), writer(std::make_shared<GraphPVDWriter>(comm, "./", "simple_linearized_solver")) {
  set_tau(tau);
}

void SimpleLinearizedSolver::set_tau(double ptau) {
  tau = ptau;

  // test values
  const double vessel_length = 1.2;
  const double radius = 0.117;
  const double wall_thickness = 0.029;
  const double elastic_modulus = 400000.0;
  const double gamma = 2;
  const double density = 1.028e-3;

  const size_t num_micro_edges = 40;

  const size_t degree = 2;

  graph = std::make_shared<GraphStorage>();

  v0 = graph->create_vertex();
  v1 = graph->create_vertex();
  v2 = graph->create_vertex();
  v3 = graph->create_vertex();
  edge1 = graph->connect(*v0, *v1, num_micro_edges);
  edge2 = graph->connect(*v1, *v2, num_micro_edges);
  edge3 = graph->connect(*v2, *v3, num_micro_edges);

  edge1->add_embedding_data({{Point(0 * vessel_length, 0, 0), Point(1 * vessel_length, 0, 0)}});
  edge2->add_embedding_data({{Point(1 * vessel_length, 0, 0), Point(2 * vessel_length, 0, 0)}});
  edge3->add_embedding_data({{Point(2 * vessel_length, 0, 0), Point(3 * vessel_length, 0, 0)}});

  auto physical_data = PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);

  edge1->add_physical_data(physical_data);
  edge2->add_physical_data(physical_data);
  edge3->add_physical_data(physical_data);

  v0->set_to_inflow_with_fixed_pressure([](double t) { return 2 * std::abs(std::sin(M_PI * t)); });
  v3->set_to_free_outflow();

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
    return get_result(*v1, *edge1);
  } else if (outlet == Outlet::out) {
    return get_result(*v2, *edge1);
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