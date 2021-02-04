////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "dof_map_network.hpp"
#include "graph_storage.hpp"

#include <cassert>

namespace macrocirculation {

SimpleDofMapNetwork::SimpleDofMapNetwork(std::size_t num_components, std::size_t num_basis_functions, std::size_t num_edges)
    : d_num_components(num_components), d_num_basis_functions(num_basis_functions), d_num_edges(num_edges) {}

void SimpleDofMapNetwork::dof_indices(const Edge &edge, std::vector<std::size_t> &dof_indices) const {
  dof_indices.resize(d_num_basis_functions * d_num_components);
  auto id = edge.get_id();
  for (auto idx = 0; idx < dof_indices.size(); idx += 1)
    dof_indices[idx] = id * d_num_basis_functions * d_num_components + idx;
}

void SimpleDofMapNetwork::dof_indices(const Edge &edge, std::vector<std::size_t> &dof_indices, std::size_t component) const {
  dof_indices.resize(d_num_basis_functions);
  auto id = edge.get_id();
  for (auto idx = 0; idx < dof_indices.size(); idx += 1)
    dof_indices[idx] = id * d_num_basis_functions * d_num_components + d_num_basis_functions*component + idx;
}

std::size_t SimpleDofMapNetwork::num_dof() const {
  return d_num_components * d_num_basis_functions * d_num_edges;
}

void extract_dof(const std::vector<std::size_t>& dof_indices, const std::vector<double>& global, std::vector<double>& local)
{
  assert(global.size() > local.size());
  assert(dof_indices.size() == local.size());

  for (std::size_t i=0; i<dof_indices.size(); i+=1)
    local[i] = global[dof_indices[i]];
}

} // namespace macrocirculation
