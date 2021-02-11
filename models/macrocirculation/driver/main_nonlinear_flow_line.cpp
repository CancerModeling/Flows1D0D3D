////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <cmath>
#include <graph_data_writer.hpp>
#include <memory>

#include "../systems/explicit_nonlinear_flow_solver.h"
#include "../systems/graph_storage.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  const double t_end = 2.;
  const std::size_t max_iter = 160000;

  const double tau = 2.5e-4/4;
  const double tau_out = 1e-3;
  const auto output_interval = static_cast< std::size_t > (tau_out/tau);

  // create the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  std::size_t ascending_aorta_id = 1;

//  auto v0 = graph->create_vertex(lm::Point(0, 0, 0));
//  auto v1 = graph->create_vertex(lm::Point(4.0, 0, 0));
//  graph->connect(*v0, *v1, ascending_aorta_id);
//  graph->refine(5);

  const std::size_t N = 11;
  auto v_prev = graph->create_vertex(lm::Point(0, 0, 0));
  v_prev->set_inflow(true);
  for (std::size_t k = 0; k < N; k += 1) {
   auto v_next = graph->create_vertex(lm::Point(4*(k + 1.) / N, 0, 0));
   graph->connect(*v_prev, *v_next, ascending_aorta_id);
   v_prev = v_next;
  }
  std::cout << graph->num_vertices() << std::endl;

  mc::ExplicitNonlinearFlowSolver solver(graph);
  solver.set_tau(tau);
  solver.use_ssp_method();

  std::vector< double > Q_vertex_values(graph->num_edges() * 2, 0);
  std::vector< double > A_vertex_values(graph->num_edges() * 2, 0);

  for (std::size_t it = 0; it < max_iter; it+=1)
  {
    std::cout << "iter " << it << std::endl;

    solver.solve();

    if (it % output_interval == 0)
    {
      // save solution
      solver.get_solution_on_vertices(Q_vertex_values, A_vertex_values);
      mc::GraphDataWriter writer;
      writer.add_vertex_data("Q", Q_vertex_values);
      writer.add_vertex_data("A", A_vertex_values);
      writer.write_vtk("line_solution", *graph, it);
    }

    // break
    if (solver.get_time() > t_end+1e-12)
      break;
  }
}
