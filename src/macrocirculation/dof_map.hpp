////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_DOF_MAP_HPP
#define TUMORMODELS_DOF_MAP_HPP

#include <functional>
#include <memory>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class Vertex;
class Edge;
class GraphStorage;
class MicroEdge;
class PetscVec;

/*! @brief Simple dof map, which orders the dofs on the macro edge by the formula
 *         dof_interval_start + edge_id * num_components * num_basis_functions + num_basis_functions * component_id + basis_function_id
 */
class LocalEdgeDofMap {
public:
  LocalEdgeDofMap(std::size_t dof_interval_start,
                  std::size_t num_components,
                  std::size_t num_basis_functions,
                  std::size_t num_micro_edges);

  void dof_indices(std::size_t local_micro_edge_id, std::size_t component, std::vector<std::size_t> &dof_indices) const;

  void dof_indices(const MicroEdge &micro_edge, std::size_t component, std::vector<std::size_t> &dof_indices) const;

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

/** @brief Simple dof map, which saves the dofs on the macro vertex. These could for instance belong to 0D model. */
class LocalVertexDofMap {
public:
  LocalVertexDofMap(std::size_t dof_interval_start, std::size_t num_components);

  std::size_t num_local_dof() const;

  const std::vector<std::size_t> &dof_indices() const;

private:
  std::vector<std::size_t> d_dof_indices;
};

/** @brief Stores and returns the local dof-maps for all macro primitives. */
class DofMap {
public:
  explicit DofMap(const GraphStorage &graph);

  DofMap(std::size_t num_vertices, std::size_t num_edges);

  const LocalEdgeDofMap &get_local_dof_map(const Edge &e) const;

  const LocalVertexDofMap &get_local_dof_map(const Vertex &e) const;

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

  /*! @brief Creates the dof-map.
   *
   * @param comm    The communicator.
   * @param graph   The graph storage
   * @param num_components   The number of components which we want to calculate.
   * @param degree           The degree of the FE-shape functions.
   * @param start_dof_offset The start offset in the dof mapping.
   * @param global           If global is true, the dofs are distributed to all macro-edges,
   *                         starting with the edges assigned to rank 0, then rank 1 and so on.
   *                         If global is false, then the dofs are only distributed to the macro-edges
   *                         which belong to the calling rank.
   *                         The first approach is useful if you want to assemble a global matrix with e.g. petsc.
   *                         The second approach is useful if you have a fully explicit scheme and only need the local dofs.
   */
  void create(MPI_Comm comm,
              const GraphStorage &graph,
              std::size_t num_components,
              std::size_t degree,
              std::size_t start_dof_offset,
              bool global);

  void create(MPI_Comm comm,
              const GraphStorage &graph,
              std::size_t num_components,
              std::size_t degree,
              std::size_t start_dof_offset,
              bool global,
              const std::function<size_t(const Vertex &)> &num_vertex_dofs);

  static void create(MPI_Comm comm,
              const std::vector< std::shared_ptr< GraphStorage > > &graphs,
              const std::vector< std::shared_ptr< DofMap > > &dof_maps,
              std::size_t num_components,
              std::size_t degree,
              const std::function<size_t(const GraphStorage&, const Vertex &)> &num_vertex_dofs);

  /*! @brief Returns the first global dof index in this dof map.
   *         If no global distribution is used the first local index (0) is returned. */
  size_t first_global_dof() const;

  /*! @brief Returns the number of dof indices in this dof map.
   *         If no global distribution is used this corresponds to the number of local dofs.
   *         If a global distribution is used this corresponds to the number of global dofs.
   */
  std::size_t num_dof() const;


  /*! @brief Returns the last global dof index in this dof map.
   *         If no global distribution is used the last local index is returned, which is equal to the number of dofs.
   */
  size_t last_global_dof() const;

  size_t first_owned_global_dof() const;

  size_t num_owned_dofs() const;

private:

  std::vector<std::unique_ptr<LocalEdgeDofMap>> d_local_dof_maps;

  std::vector<std::unique_ptr<LocalVertexDofMap>> d_local_vertex_dof_maps;

  std::size_t d_first_global_dof;
  std::size_t d_num_dof;

  size_t d_first_owned_global_dof;
  size_t d_num_owned_dofs;

private:
  void create(MPI_Comm comm, const std::vector<GraphStorage> &graphs, const std::vector<DofMap> &dof_maps, size_t num_components, size_t degree, bool global, const std::function<size_t(const Vertex &)> &num_vertex_dofs);

  void add_local_dof_map(const Edge &e,
                         std::size_t num_components,
                         std::size_t num_basis_functions,
                         std::size_t num_micro_edges,
                         std::size_t start_dof);

  void add_local_dof_map(const Vertex &v, std::size_t start_dof, std::size_t num_components);
};

/*! @brief Copies the dof values in dof_indices from a global vector into a local vector. */
void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const std::vector<double> &global,
                 std::vector<double> &local);

/*! @brief Copies the dof values in dof_indices from a global vector into a local vector. */
void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const PetscVec &global,
                 std::vector<double> &local);

} // namespace macrocirculation

#endif