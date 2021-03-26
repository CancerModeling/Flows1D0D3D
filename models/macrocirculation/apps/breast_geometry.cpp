////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <macrocirculation/dof_map.hpp>
#include <macrocirculation/explicit_nonlinear_flow_solver.hpp>
#include <macrocirculation/graph_pvd_writer.hpp>
#include <macrocirculation/interpolate_to_vertices.hpp>
#include <memory>

#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append("coarse-network-geometry.json", *graph);

  auto dof_map = std::make_shared<mc::DofMap>(graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  std::vector< double > u(dof_map->num_dof(), 0);

  std::vector< mc::Point > points;
  std::vector< double > u_vertex_values;

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "breast_geometry_solution");
  mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 0, u, points, u_vertex_values);

  pvd_writer.set_points(points);
  pvd_writer.add_vertex_data("u", u_vertex_values);
  pvd_writer.write(0);

  MPI_Finalize();
}
