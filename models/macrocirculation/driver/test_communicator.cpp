////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <mpi.h>
#include <memory>
#include <gmm.h>

#include "../systems/graph_storage.hpp"
#include "../systems/graph_partitioner.hpp"
#include "../systems/communication/mpi.hpp"
#include "../systems/dof_map_network.hpp"
#include "../systems/communicator.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  MPI_Init( &argc, &argv );

  const std::size_t num_edges_per_segment = 8;

  // create the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto start = graph->create_vertex(lm::Point(0, 0, 0));
  auto end = graph->create_vertex(lm::Point(4, 0, 0));
  graph->line_to(*start, *end, 0, num_edges_per_segment);

  // partition graph
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map = std::make_shared<mc::SimpleDofMapNetwork>(2, degree + 1, graph->num_edges());

  std::vector<double> vec_u_src(dof_map->num_dof(), 0);

  const auto rank = mc::mpi::rank(MPI_COMM_WORLD);
  const auto size = mc::mpi::size(MPI_COMM_WORLD);
  std::cout << "rank = " << rank << " size = " << size << std::endl;

  std::vector< std::size_t > dof_indices(dof_map->num_local_dof());
  for (auto edge_id: graph->get_active_edge_ids(rank))
  {
    std::cout << "rank " << rank << " edge_id " << edge_id << std::endl;
    const auto edge = graph->get_edge(edge_id);
    dof_map->dof_indices(*edge, dof_indices);
    for (auto idx: dof_indices)
      vec_u_src[idx] = idx;
  }

  std::cout << " rank = " << rank << " " << vec_u_src << std::endl;

  mc::Communicator communicator(MPI_COMM_WORLD, graph, dof_map);
  communicator.update_ghost_layer(vec_u_src);

  std::cout << " rank = " << rank << " " << vec_u_src << std::endl;

  communicator.gather(1, vec_u_src);

  std::cout << " rank = " << rank << " " << vec_u_src << std::endl;

  MPI_Finalize();
}
