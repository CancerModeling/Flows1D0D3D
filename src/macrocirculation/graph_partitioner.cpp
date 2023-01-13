////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_partitioner.hpp"
#include "communication/mpi.hpp"
#include "graph_storage.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <queue>

namespace macrocirculation {

void naive_mesh_partitioner(GraphStorage &graph, MPI_Comm comm) {
  int size;
  MPI_Comm_size(comm, &size);

  for (std::size_t rank = 0; rank < static_cast< size_t > ( size ); rank += 1) {
    std::size_t start_edge_id = (rank * graph.num_edges()) / size;
    std::size_t end_edge_id = ((rank + 1) * graph.num_edges()) / size;

    for (std::size_t edge_id = start_edge_id; edge_id < end_edge_id; edge_id += 1) {
      auto edge = graph.get_edge(edge_id);
      graph.assign_edge_to_rank(*edge, rank);
    }
  }
}

struct RankPriorityPair {
  int rank;
  int total_priority;
};

struct PrioritizedEdge {
  size_t edge_id;
  int total_priority;
};

void priority_mesh_partitioner(MPI_Comm comm, GraphStorage &graph, const std::function<int(const Edge &)> &estimator) {
  auto cmp = [](const RankPriorityPair &left, const RankPriorityPair &right) {
    return left.total_priority > right.total_priority;
  };

  std::priority_queue<RankPriorityPair, std::vector<RankPriorityPair>, decltype(cmp)> queue(cmp);

  for (auto rank = 0; rank < mpi::size(comm); rank += 1) {
    queue.push({rank, 0});
  }

  std::vector<PrioritizedEdge> edges;

  for (auto e_id : graph.get_edge_ids()) {
    auto &edge = *graph.get_edge(e_id);
    auto priority = estimator(edge);
    edges.push_back({e_id, priority});
  }

  // largest edges first
  std::sort(edges.begin(), edges.end(), [](const auto &l, const auto &r) {
    return l.total_priority > r.total_priority;
  });

  for (auto pedge : edges) {
    RankPriorityPair rank_priority = queue.top();
    queue.pop();

    auto &edge = *graph.get_edge(pedge.edge_id);
    graph.assign_edge_to_rank(edge, rank_priority.rank);

    rank_priority.total_priority += pedge.total_priority;
    queue.push(rank_priority);
  }
}

int flow_cost_esimator(const GraphStorage &graph, const Edge &e, size_t degree) {
  size_t num_dofs = e.num_micro_edges() * (degree + 1);

  for (auto v_id : e.get_vertex_neighbors()) {
    auto &vertex = *graph.get_vertex(v_id);
    if (!vertex.bc_finalized())
      throw std::runtime_error("boundary conditions need to be set before distributing the graph");
    if (vertex.is_windkessel_outflow())
      num_dofs += 1;
    if (vertex.is_vessel_tree_outflow())
      num_dofs += vertex.get_vessel_tree_data().resistances.size();
  }

  return static_cast<int>(num_dofs);
}

void flow_mesh_partitioner(MPI_Comm comm, GraphStorage &graph, size_t degree) {

  auto estimator = [&](const Edge &e) {
    return flow_cost_esimator(graph, e, degree);
  };
  priority_mesh_partitioner(comm, graph, estimator);
}

} // namespace macrocirculation
