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

namespace macrocirculation {

LocalDofMap::LocalDofMap(std::size_t dof_interval_start,
                         std::size_t num_components,
                         std::size_t num_basis_functions,
                         std::size_t num_micro_edges)
    : d_dof_interval_start(dof_interval_start),
      d_dof_interval_end(dof_interval_start + num_micro_edges * num_components * num_basis_functions),
      d_num_components(num_components),
      d_num_basis_functions(num_basis_functions),
      d_num_micro_edges(num_micro_edges) {}

void LocalDofMap::dof_indices(std::size_t micro_edge_index,
                              std::size_t component,
                              std::vector<std::size_t> &dof_indices) const {
  assert(dof_indices.size() == d_num_basis_functions);
  for (auto idx = 0; idx < d_num_basis_functions; idx += 1)
    dof_indices[idx] = d_dof_interval_start + micro_edge_index * d_num_basis_functions * d_num_components +
                       d_num_basis_functions * component + idx;
  assert(dof_indices.back() < d_dof_interval_end);
}

std::size_t LocalDofMap::num_local_dof() const {
  return d_num_components * d_num_basis_functions * d_num_micro_edges;
}

std::size_t LocalDofMap::num_micro_edges() const {
  return d_num_micro_edges;
}

std::size_t LocalDofMap::num_micro_vertices() const {
  return d_num_micro_edges + 1;
}

std::size_t LocalDofMap::num_components() const {
  return d_num_components;
}

std::size_t LocalDofMap::num_basis_functions() const {
  return d_num_basis_functions;
}

DofMap::DofMap(std::size_t num_edges)
    : d_local_dof_maps(num_edges), d_num_dof(0) {}

void DofMap::add_local_dof_map(const Edge &e,
                               std::size_t num_components,
                               std::size_t num_basis_functions,
                               std::size_t num_local_micro_edges) {
  assert(d_local_dof_maps[e.get_id()] == nullptr);
  auto local_dof_map = std::make_unique<LocalDofMap>(d_num_dof, num_components, num_basis_functions, num_local_micro_edges);
  d_num_dof += local_dof_map->num_local_dof();
  d_local_dof_maps[e.get_id()] = std::move(local_dof_map);
}

const LocalDofMap &DofMap::get_local_dof_map(const Edge &e) const {
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

void DofMap::create_for_node(MPI_Comm comm,
                    const GraphStorage &graph,
                    std::size_t num_components,
                    std::size_t degree) {
  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);
    add_local_dof_map(*edge, num_components, degree + 1, edge->num_micro_edges());
  }
}

} // namespace macrocirculation
