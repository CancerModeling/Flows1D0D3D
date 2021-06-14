////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_partitioner.hpp"
#include "graph_storage.hpp"
#include <cmath>
#include <mpi.h>

namespace macrocirculation {

void naive_mesh_partitioner(GraphStorage &graph, MPI_Comm comm) {
  int size;
  MPI_Comm_size(comm, &size);

  for (std::size_t rank = 0; rank < size; rank += 1) {
    std::size_t start_edge_id = (rank * graph.num_edges()) / size;
    std::size_t end_edge_id = ((rank + 1) * graph.num_edges()) / size;

    for (std::size_t edge_id = start_edge_id; edge_id < end_edge_id; edge_id += 1) {
      auto edge = graph.get_edge(edge_id);
      graph.assign_edge_to_rank(*edge, rank);
    }
  }
}

} // namespace macrocirculation
