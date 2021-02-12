////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <cmath>
#include <graph_data_writer.hpp>
#include <memory>

#include "../systems/explicit_nonlinear_flow_solver.hpp"
#include "../systems/graph_storage.hpp"
#include "../systems/vessel_data_storage.hpp"
#include "../systems/vessel_formulas.hpp"

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

  const std::size_t num_edges_per_segment = 11;

  // we create data for the ascending aorta
  auto vessel_data = std::make_shared< mc::VesselDataStorage > ();
  std::size_t ascending_aorta_id = vessel_data->add_parameter({
    592.4e2, // 592.4 10^2 Pa,  TODO: Check if units are consistent!
    6.97,    // 6.97 cm^2,      TODO: Check if units are consistent!
    1.028,   // 1.028 kg/cm^3,  TODO: Check if units are consistent!
  });

  // create the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto start = graph->create_vertex(lm::Point(0, 0, 0));
  graph->line_to(start, lm::Point(4, 0, 0), ascending_aorta_id, num_edges_per_segment);

  // set inflow boundary conditions
  start->set_to_inflow(mc::heart_beat_inflow);

  // configure solver
  mc::ExplicitNonlinearFlowSolver solver(graph, vessel_data);
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
