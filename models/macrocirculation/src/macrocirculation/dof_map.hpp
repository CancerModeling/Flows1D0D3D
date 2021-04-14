////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_DOF_MAP_HPP
#define TUMORMODELS_DOF_MAP_HPP

#include <memory>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class Vertex;
class Edge;
class GraphStorage;
class MicroEdge;

/*! @brief Simple dof map, which orders the dofs by the formula
 *         dof_interval_start + edge_id * num_components * num_basis_functions + num_basis_functions * component_id + basis_function_id
 */
class LocalEdgeDofMap {
public:
  LocalEdgeDofMap(std::size_t dof_interval_start,
              std::size_t num_components,
              std::size_t num_basis_functions,
              std::size_t num_micro_edges);

  void dof_indices(std::size_t local_micro_edge_id, std::size_t component, std::vector<std::size_t> &dof_indices) const;

  void dof_indices(const MicroEdge& micro_edge, std::size_t component, std::vector<std::size_t> &dof_indices) const;

  std::size_t num_local_dof() const;
  std::size_t num_micro_edges() const;
  std::size_t num_micro_vertices() const;
  std::size_t num_components() const;
  std::size_t num_basis_functions() const;

private:
  std::size_t d_dof_interval_start;
  std::size_t d_dof_interval_end;

  std::size_t d_num_components;
  std::size_t d_num_basis_functions;
  std::size_t d_num_micro_edges;
};

class LocalVertexDofMap {
public:
  LocalVertexDofMap(std::size_t dof_interval_start, std::size_t num_components);

  std::size_t num_local_dof() const;

  const std::vector< std::size_t > & dof_indices() const;

private:
  std::vector< std::size_t > d_dof_indices;
};

class DofMap {
public:
  explicit DofMap(std::size_t num_vertices, std::size_t num_edges);

  const LocalEdgeDofMap &get_local_dof_map(const Edge &e) const;

  const LocalVertexDofMap &get_local_dof_map(const Vertex &e) const;

  void add_local_dof_map(const Edge &e,
                         std::size_t num_components,
                         std::size_t num_basis_functions,
                         std::size_t num_micro_edges);

  void add_local_dof_map(const Vertex &v, std::size_t num_components);

  std::size_t num_dof() const;

  /*! @brief Creates the dof-map.
   *
   * @param comm    The communicator.
   * @param graph   The graph storage
   * @param num_components  The number of components which we want to calculate.
   * @param degree          The degree of the FE-shape functions.
   * @param global          If global is true, the dofs are distributed to all macro-edges,
   *                        starting with the edges assigned to rank 0, then rank 1 and so on.
   *                        If global is false, then the dofs are only distributed to the macro-edges
   *                        which belong to the calling rank.
   *                        The first approach is useful if you want to assemble a global matrix with e.g. petsc.
   *                        The second approach is useful if you have a fully explicit scheme and only need the local dofs.
   */
  void create(MPI_Comm comm,
              const GraphStorage &graph,
              std::size_t num_components,
              std::size_t degree,
              bool global);

private:
  std::vector<std::unique_ptr<LocalEdgeDofMap>> d_local_dof_maps;

  std::vector<std::unique_ptr<LocalVertexDofMap>> d_local_vertex_dof_maps;

  std::size_t d_num_dof;
};

/*! @brief Copies the dof values in dof_indices from a global vector into a local vector. */
void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const std::vector<double> &global,
                 std::vector<double> &local);

} // namespace macrocirculation

#endif