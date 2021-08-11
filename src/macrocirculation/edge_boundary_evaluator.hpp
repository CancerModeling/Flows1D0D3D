////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EDGE_BOUNDARY_EVALUATOR_HPP
#define TUMORMODELS_EDGE_BOUNDARY_EVALUATOR_HPP

#include <memory>
#include <mpi.h>
#include <vector>

#include "communicator.hpp"

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class Vertex;
class Edge;
class PetscVec;

/*!@brief Evaluates a finite element function at the macro edge boundaries
 *        and communicates the values to all ranks with an adjacent primitive.
 */
class EdgeBoundaryEvaluator {
public:
  EdgeBoundaryEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t component);

  /*! @brief Evaluates and communicates the values of at the macro edge boundaries. */
  void init(const std::vector<double> &u_prev);

  /*! @brief Evaluates and communicates the values of at the macro edge boundaries. */
  void init(const PetscVec &u_prev);

  /*! @brief Returns the boundary value at the given vertex on the given macro edge. */
  double operator()(const Vertex &vertex, const Edge &edge) const;

  /*! @brief Returns all the boundary values at the given vertex for all macro edge neighbors. */
  void operator()(const Vertex &vertex, std::vector<double> &values) const;

  /*! @brief Returns the two boundary values at the given macro edge. */
  void operator()(const Edge &edge, std::vector<double> &values) const;

private:
  template<typename VectorType>
  void evaluate_macro_edge_boundary_values(const VectorType &u_prev);

private:
  MPI_Comm d_comm;

  /*! @brief The current domain on which we evaluate the given function. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain. */
  std::shared_ptr<DofMap> d_dof_map;

  /*! @brief The dof component we want to evaluate at the boundaries. */
  size_t d_component;

  /*! @brief Communicates the evaluated values to other ranks. */
  Communicator d_edge_boundary_communicator;

  /*! @brief Contains the values of the given component at the left and right boundary point of the given macro-edge.
    *         For the 2i-th edge the left boundary value is at 2*i, while the right entry is at 2*i+1.
    */
  std::vector<double> d_macro_edge_boundary_value;
};

} // namespace macrocirculation

#endif