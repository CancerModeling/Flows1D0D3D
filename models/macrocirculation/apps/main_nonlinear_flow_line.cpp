////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <memory>

#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_data_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

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
  auto end = graph->create_vertex(lm::Point(4, 0, 0));
  graph->line_to(*start, *end, ascending_aorta_id, num_edges_per_segment);

  // set inflow boundary conditions
  start->set_to_inflow(mc::heart_beat_inflow());

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> solver(MPI_COMM_WORLD, graph, vessel_data);
  solver.set_tau(tau);
  solver.use_ssp_method();

  std::vector< double > Q_vertex_values(graph->num_edges() * 2, 0);
  std::vector< double > A_vertex_values(graph->num_edges() * 2, 0);
  std::vector< double > p_static_vertex_values(graph->num_edges() * 2, 0);
  std::vector< double > p_total_vertex_values(graph->num_edges() * 2, 0);
  std::vector< libMesh::Point > points(graph->num_edges() * 2);

  mc::GraphPVDWriter writer(MPI_COMM_WORLD, "./", "line_solution");

  for (std::size_t it = 0; it < max_iter; it+=1)
  {
    std::cout << "iter " << it << std::endl;

    solver.solve();

    if (it % output_interval == 0)
    {
      // save solution
      solver.get_solution_on_vertices(points, Q_vertex_values, A_vertex_values);
      solver.get_total_pressure_on_vertices(points, p_total_vertex_values);
      solver.get_static_pressure_on_vertices(points, p_static_vertex_values);
      writer.set_points(points);
      writer.add_vertex_data("Q", Q_vertex_values);
      writer.add_vertex_data("A", A_vertex_values);
      writer.add_vertex_data("p_total", p_total_vertex_values);
      writer.add_vertex_data("p_static", p_static_vertex_values);
      writer.write(solver.get_time());
    }

    // break
    if (solver.get_time() > t_end+1e-12)
      break;
  }
}
