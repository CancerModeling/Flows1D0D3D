////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <memory>
#include "libmesh/libmesh.h"
#include <cmath>

#include "macrocirculation/implicit_advection_solver.hpp"
#include "macrocirculation/graph_storage.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 100;

  const double velocity = 1;
  const double tau = 0.001;
  const double t_end = 2;
  const auto inflow_boundary_value = [](double t_now) -> double { return std::sin(M_PI * 3 * t_now); };
  // const auto inflow_boundary_value = [](double t_now) -> double { return 1; };

  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  // create the ascending aorta
  auto graph = std::make_shared< mc::GraphStorage >();

  std::size_t ascending_aorta_id = 1;

  auto v0 = graph->create_vertex();
  auto v1 = graph->create_vertex();
  auto edge = graph->connect(*v0, *v1, num_micro_edges);
  edge->add_embedding_data({{ mc::Point(0, 0, 0), mc::Point(1, 0, 0) }});
  edge->add_physical_data({ 0, 0, 0, 1. });

  // const std::size_t N = 32;
  // auto v_prev = graph.create_vertex(lm::Point(0, 0, 0));
  // for (std::size_t k = 0; k < N; k += 1) {
  //  auto v_next = graph.create_vertex(lm::Point((k + 1.) / N, 0, 0));
  //  graph.connect(*v_prev, *v_next, ascending_aorta_id);
  //  v_prev = v_next;
  // }

  mc::ImplicitAdvectionSolver solver(graph, degree, tau, t_end, velocity, inflow_boundary_value);
  solver.set_output_interval(100);

  solver.solve();
}