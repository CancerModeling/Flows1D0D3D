////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_COMMUNICATOR_HPP
#define TUMORMODELS_COMMUNICATOR_HPP

#include <memory>
#include <mpi.h>
#include <vector>

#include "communication/buffer.hpp"

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMapNetwork;

class Communicator {
public:
  Communicator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMapNetwork> dof_map);

  /*! @brief Updates the ghost layer on all ranks for the given vector. */
  void update_ghost_layer(std::vector<double> &u);

  /*! @brief Gathers all the dof data on one single rank. */
  void gather(int rank, std::vector<double> &u);

private:
  /*! @brief Unpacks the dof already in the receive buffer into the given vector. */
  void unpack_dof_from_receive_buffer(std::vector<double> &u);

  /*! @brief Packs the dof on the edges with the given ids into the send buffer of the receiving rank. */
  void pack_dof_into_send_buffer(int receiver_rank, const std::vector<std::size_t> &edge_ids, const std::vector<double> &u);

  MPI_Comm d_comm;
  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMapNetwork> d_dof_map;
  BufferSystem d_buffer_system;
  int d_rank;
  std::vector<std::vector<std::size_t>> d_ghost_dofs_to_send;
};

} // namespace macrocirculation


#endif //TUMORMODELS_COMMUNICATOR_HPP
