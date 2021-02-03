////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <memory>
#include "libmesh/libmesh.h"
#include <cmath>

#include "../systems/advection_solver.hpp"
#include "../systems/graph_storage.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  constexpr std::size_t degree = 3;

  const double velocity = 2;
  const double tau = 0.1;
  const double t_end = 1;
  const auto inflow_boundary_value = [](double t_now) -> double { return std::sin(M_PI * t_now); };

  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  // create the ascending aorta
  auto graph = std::make_shared< mc::GraphStorage >();

  std::size_t ascending_aorta_id = 1;

  auto v0 = graph->create_vertex(lm::Point(0, 0, 0));
  auto v1 = graph->create_vertex(lm::Point(1.0, 0, 0));
  graph->connect(*v0, *v1, ascending_aorta_id);
  graph->refine(5);

  // const std::size_t N = 32;
  // auto v_prev = graph.create_vertex(lm::Point(0, 0, 0));
  // for (std::size_t k = 0; k < N; k += 1) {
  //  auto v_next = graph.create_vertex(lm::Point((k + 1.) / N, 0, 0));
  //  graph.connect(*v_prev, *v_next, ascending_aorta_id);
  //  v_prev = v_next;
  // }

  mc::AdvectionSolver<degree> solver(graph, tau, t_end, velocity, inflow_boundary_value);
}
