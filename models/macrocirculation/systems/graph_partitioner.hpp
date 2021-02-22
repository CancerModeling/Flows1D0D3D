////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_PARTITIONER_HPP
#define TUMORMODELS_GRAPH_PARTITIONER_HPP

#include <mpi.h>

namespace macromesh {

// forward declarations
class GraphStorage;

/*! @brief Distributes a graph to several processors through the ids of the assigned primitives.
 *         This partitioning makes sense when the edge ids and the geometric position correlate.
 */
void naive_mesh_partitioner(GraphStorage & graph, MPI_Comm comm=MPI_COMM_WORLD);

}

#endif //TUMORMODELS_GRAPH_PARTITIONER_HPP
