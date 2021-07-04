////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_PARTITIONER_HPP
#define TUMORMODELS_GRAPH_PARTITIONER_HPP

#include <functional>
#include <mpi.h>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class Edge;

/*! @brief Distributes a graph to several processors through the ids of the assigned primitives.
 *         This partitioning makes sense when the edge ids and the geometric position correlate.
 */
void naive_mesh_partitioner(GraphStorage &graph, MPI_Comm comm = MPI_COMM_WORLD);

/*! @brief Distributes a graph to several processors using a cost estimator function for the primitives.
 */
void priority_mesh_partitioner(MPI_Comm comm, GraphStorage &graph, const std::function<int (const Edge&)> & estimator);

/*! @brief Distributes a graph to several processors using a cost estimator for the flow problem;
 */
void flow_mesh_partitioner(MPI_Comm comm, GraphStorage &graph, size_t degree);

} // namespace macrocirculation

#endif //TUMORMODELS_GRAPH_PARTITIONER_HPP
