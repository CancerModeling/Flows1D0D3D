////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include "mpi.h"
#include <cmath>
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"

namespace mc = macrocirculation;

/*! @brief Checks if we converge to a constant solution if we impose fixed characteristic boundary conditions. */
TEST_CASE("LinearCharacteristicBCs", "[LinearCharacteristicBCs]") {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 20;

  const double tau = 1e-3;
  const double t_end = 0.17;

  // vessel parameters for a given vessel
  const double vessel_length = 42.2;
  const double radius = 0.403;
  const double wall_thickness = 0.067;
  const double elastic_modulus = 400000.0;
  const double gamma = 9;
  const double density = 1.028e-3;

  auto graph = std::make_shared<mc::GraphStorage>();

  auto v0 = graph->create_vertex();
  auto v1 = graph->create_vertex();
  auto edge1 = graph->connect(*v0, *v1, num_micro_edges);

  auto physical_data = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
  physical_data.viscosity = 0.;

  edge1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(2, 0, 0)}});
  edge1->add_physical_data(physical_data);

  const double p_in = 5.;
  const double q_in = 4.;

  v0->set_to_linear_characteristic_inflow(mc::LinearFlowSolver::get_C(*edge1), mc::LinearFlowSolver::get_L(*edge1), true, p_in, q_in);
  v1->set_to_linear_characteristic_inflow(mc::LinearFlowSolver::get_C(*edge1), mc::LinearFlowSolver::get_L(*edge1), false, p_in, q_in);

  mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(PETSC_COMM_WORLD, *graph, 2, degree, true);

  mc::LinearFlowSolver solver(PETSC_COMM_WORLD, graph, dof_map, degree);
  solver.setup(tau);

  double t = 0;
  const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));

  for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
    t += tau;
    solver.solve(tau, t);
  }

  double p, q;
  // at the left tip
  solver.evaluate_1d_pq_values(*edge1, 0, p, q);
  REQUIRE(p == Approx(p_in).margin(3e-4));
  REQUIRE(q == Approx(q_in).margin(3e-4));

  // in the middle
  solver.evaluate_1d_pq_values(*edge1, 0.5, p, q);
  REQUIRE(p == Approx(p_in).margin(3e-4));
  REQUIRE(q == Approx(q_in).margin(3e-4));

  // at the right tip
  solver.evaluate_1d_pq_values(*edge1, 1, p, q);
  REQUIRE(p == Approx(p_in).margin(3e-4));
  REQUIRE(q == Approx(q_in).margin(3e-4));
}
