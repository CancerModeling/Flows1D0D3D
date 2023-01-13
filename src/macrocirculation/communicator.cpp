////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "communicator.hpp"

#include "communication/mpi.hpp"
#include "graph_storage.hpp"
#include <utility>

namespace macrocirculation {

Communicator::Communicator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, DoFFunctional dof_to_send, DoFFunctional dof_to_receive)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_to_send(std::move(dof_to_send)),
      d_dof_to_receive(std::move(dof_to_receive)),
      d_buffer_system(comm, 0),
      d_rank(mpi::rank(comm)),
      d_ghost_dofs_to_send(mpi::size(comm), std::vector<std::size_t>()) {
  for (int other_rank = 0; other_rank < mpi::size(comm); other_rank += 1) {
    // we dont send to ourselves.
    if (other_rank == d_rank)
      continue;

    d_ghost_dofs_to_send[other_rank] = d_graph->get_ghost_edge_ids(other_rank, d_rank);
  }
}

void Communicator::update_ghost_layer(std::vector<double> &u) {
  d_buffer_system.clear();

  // send the ghost layer to our neighbors
  for (int other_rank = 0; other_rank < static_cast< int >( d_ghost_dofs_to_send.size() ); other_rank += 1) {
    // we dont send to ourselves.
    if (other_rank == d_rank)
      continue;

    pack_dof_into_send_buffer(other_rank, d_ghost_dofs_to_send[other_rank], u);
  }
  d_buffer_system.start_communication();
  d_buffer_system.end_communication();

  // receive the ghost layer from our neighbors
  unpack_dof_from_receive_buffer(u);

  // clear the buffer system to check if any information was lost
  d_buffer_system.clear();
}

void Communicator::unpack_dof_from_receive_buffer(std::vector<double> &u) {
  // receive the ghost layer from our neighbors
  for (int other_rank = 0; other_rank < static_cast< int >( d_ghost_dofs_to_send.size() ); other_rank += 1) {
    // we dont receive from ourselves.
    if (other_rank == d_rank)
      continue;

    auto &receive_buffer = d_buffer_system.get_receive_buffer(other_rank);

    while (!receive_buffer.empty()) {
      std::size_t edge_id;
      receive_buffer >> edge_id;
      assert(edge_id < d_graph->num_edges());
      const auto edge = d_graph->get_edge(edge_id);
      assert(edge->rank() == other_rank);
      const auto dof_indices = d_dof_to_receive(*edge);
      for (std::size_t idx : dof_indices)
        receive_buffer >> u[idx];
    }
  }
}

void Communicator::pack_dof_into_send_buffer(int receiver_rank, const std::vector<std::size_t> &edge_ids, const std::vector<double> &u) {
  assert(receiver_rank != d_rank);

  auto &send_buffer = d_buffer_system.get_send_buffer(receiver_rank);

  for (const auto &edge_id : edge_ids) {
    const auto edge = d_graph->get_edge(edge_id);
    const auto dof_indices = d_dof_to_send(*edge);

    send_buffer << edge_id;
    for (std::size_t idx : dof_indices)
      send_buffer << u[idx];
  }
}

Communicator Communicator::create_edge_boundary_value_communicator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph) {
  DoFFunctional boundary_edge_dofs = [](const Edge &edge) -> std::vector<std::size_t> {
    return {
      2 * edge.get_id(),    // index to the left boundary value of the given macro-edge
      2 * edge.get_id() + 1 // index to the right boundary value of the given boundary edge
    };
  };

  return Communicator{comm, std::move(graph), boundary_edge_dofs, boundary_edge_dofs};
}

} // namespace macrocirculation
