////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "dof_map.hpp"

#include <cassert>

#include "communication/mpi.hpp"
#include "graph_storage.hpp"
#include "petsc/petsc_vec.hpp"

namespace macrocirculation {

LocalEdgeDofMap::LocalEdgeDofMap(std::size_t dof_interval_start,
                                 std::size_t num_components,
                                 std::size_t num_basis_functions,
                                 std::size_t num_micro_edges)
    : d_dof_interval_start(dof_interval_start),
      d_dof_interval_end(dof_interval_start + num_micro_edges * num_components * num_basis_functions),
      d_num_components(num_components),
      d_num_basis_functions(num_basis_functions),
      d_num_micro_edges(num_micro_edges) {}

void LocalEdgeDofMap::dof_indices(std::size_t micro_edge_index,
                                  std::size_t component,
                                  std::vector<std::size_t> &dof_indices) const {
  assert(dof_indices.size() == d_num_basis_functions);
  for (auto idx = 0; idx < d_num_basis_functions; idx += 1)
    dof_indices[idx] = d_dof_interval_start + micro_edge_index * d_num_basis_functions * d_num_components +
                       d_num_basis_functions * component + idx;
  assert(dof_indices.back() < d_dof_interval_end);
}

void LocalEdgeDofMap::dof_indices(const MicroEdge &micro_edge, std::size_t component, std::vector<std::size_t> &dof_indices_vec) const {
  dof_indices(micro_edge.get_local_id(), component, dof_indices_vec);
}

std::size_t LocalEdgeDofMap::num_local_dof() const {
  return d_num_components * d_num_basis_functions * d_num_micro_edges;
}

std::size_t LocalEdgeDofMap::num_micro_edges() const {
  return d_num_micro_edges;
}

std::size_t LocalEdgeDofMap::num_micro_vertices() const {
  return d_num_micro_edges + 1;
}

std::size_t LocalEdgeDofMap::num_components() const {
  return d_num_components;
}

std::size_t LocalEdgeDofMap::num_basis_functions() const {
  return d_num_basis_functions;
}

DofMap::DofMap(std::size_t num_vertices, std::size_t num_edges)
    : d_local_dof_maps(num_edges),
      d_local_vertex_dof_maps(num_vertices),
      d_num_dof(0),
      d_first_owned_global_dof(0),
      d_num_owned_dofs(0) {}

void DofMap::add_local_dof_map(const Edge &e,
                               std::size_t num_components,
                               std::size_t num_basis_functions,
                               std::size_t num_local_micro_edges) {
  assert(d_local_dof_maps[e.get_id()] == nullptr);
  auto local_dof_map = std::make_unique<LocalEdgeDofMap>(d_num_dof, num_components, num_basis_functions, num_local_micro_edges);
  d_num_dof += local_dof_map->num_local_dof();
  d_local_dof_maps[e.get_id()] = std::move(local_dof_map);
}

void DofMap::add_local_dof_map(const Vertex &v, std::size_t num_components) {
  assert(d_local_vertex_dof_maps[v.get_id()] == nullptr);
  auto local_dof_map = std::make_unique<LocalVertexDofMap>(d_num_dof, num_components);
  d_num_dof += local_dof_map->num_local_dof();
  d_local_vertex_dof_maps[v.get_id()] = std::move(local_dof_map);
}

const LocalEdgeDofMap &DofMap::get_local_dof_map(const Edge &e) const {
  if (d_local_dof_maps.at(e.get_id()) == nullptr)
    throw std::runtime_error("dof map for edge with id " + std::to_string(e.get_id()) + " was not initialized.");
  return *d_local_dof_maps.at(e.get_id());
}

std::size_t DofMap::num_dof() const {
  return d_num_dof;
}

void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const std::vector<double> &global,
                 std::vector<double> &local) {
  assert(global.size() > local.size());
  assert(dof_indices.size() == local.size());

  for (std::size_t i = 0; i < dof_indices.size(); i += 1)
    local[i] = global[dof_indices[i]];
}

void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const PetscVec &global,
                 std::vector<double> &local) {
  assert(dof_indices.size() == local.size());

  for (std::size_t i = 0; i < dof_indices.size(); i += 1)
    local[i] = global.get(dof_indices[i]);
}

void DofMap::create(MPI_Comm comm,
                    const GraphStorage &graph,
                    std::size_t num_components,
                    std::size_t degree,
                    bool global) {
  for (int rank = 0; rank < mpi::size(comm); rank += 1) {
    const bool callingRank = (rank == mpi::rank(comm));

    // if we do not assign all the ranks AND this is not the calling rank we add nothing
    if (!global && !callingRank)
      continue;

    // if we are the rank
    if (callingRank)
      d_first_owned_global_dof = d_num_dof;

    for (const auto &e_id : graph.get_active_edge_ids(rank)) {
      const auto edge = graph.get_edge(e_id);
      add_local_dof_map(*edge, num_components, degree + 1, edge->num_micro_edges());
      if (callingRank)
        d_num_owned_dofs += get_local_dof_map(*edge).num_local_dof();
    }

    for (const auto &v_id : graph.get_active_vertex_ids(rank)) {
      const auto vertex = graph.get_vertex(v_id);

      if (!vertex->bc_finalized())
        throw std::runtime_error(
          "Boundary conditions have to be finalized before distributing the dof on primitives.\n"
          "Please call GraphStorage::finalize_bcs() before.");


      if (vertex->is_windkessel_outflow()) {
        add_local_dof_map(*vertex, 1);
        if (callingRank)
          d_num_owned_dofs += get_local_dof_map(*vertex).num_local_dof();
      } else if (vertex->is_vessel_tree_outflow()) {
        add_local_dof_map(*vertex, vertex->get_vessel_tree_data().resistances.size());
        if (callingRank)
          d_num_owned_dofs += get_local_dof_map(*vertex).num_local_dof();
      } else if (vertex->is_rcl_outflow()) {
        add_local_dof_map(*vertex, 2 * vertex->get_rcl_data().resistances.size());
        if (callingRank)
          d_num_owned_dofs += get_local_dof_map(*vertex).num_local_dof();
      }
    }
  }
}

const LocalVertexDofMap &DofMap::get_local_dof_map(const Vertex &v) const {
  if (d_local_vertex_dof_maps.at(v.get_id()) == nullptr)
    throw std::runtime_error("dof map for vertex with id " + std::to_string(v.get_id()) + " was not initialized.");
  return *d_local_vertex_dof_maps.at(v.get_id());
}

LocalVertexDofMap::LocalVertexDofMap(std::size_t dof_interval_start, std::size_t num_components) {
  for (size_t c = 0; c < num_components; c += 1)
    d_dof_indices.push_back(dof_interval_start + c);
}

size_t DofMap::first_owned_global_dof() const { return d_first_owned_global_dof; }

size_t DofMap::num_owned_dofs() const { return d_num_owned_dofs; }

std::size_t LocalVertexDofMap::num_local_dof() const { return d_dof_indices.size(); }

const std::vector<std::size_t> &LocalVertexDofMap::dof_indices() const { return d_dof_indices; }

} // namespace macrocirculation
