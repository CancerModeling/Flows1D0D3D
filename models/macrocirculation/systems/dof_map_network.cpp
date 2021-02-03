////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "dof_map_network.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

SimpleDofMapNetwork::SimpleDofMapNetwork(std::size_t num_components, std::size_t num_basis_functions)
    : d_num_components(num_components), d_num_basis_functions(num_basis_functions) {}

void SimpleDofMapNetwork::dof_indices(const Edge &edge, std::vector<std::size_t> &dof_indices) const {
  dof_indices.resize(d_num_basis_functions * d_num_components);
  auto id = edge.get_id();
  for (auto idx = 0; idx < dof_indices.size(); idx += 1)
    dof_indices[idx] = id * d_num_basis_functions * d_num_components + idx;
}

} // namespace macrocirculation
