////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_DOF_MAP_NETWORK_HPP
#define TUMORMODELS_DOF_MAP_NETWORK_HPP

#include <vector>

namespace macrocirculation {

// forward declarations
class Edge;

/*! @brief Abstract interface for the dof map in the 1D network. */
class DofMapNetwork {
public:
  virtual ~DofMapNetwork() = default;

  virtual void dof_indices(const Edge &edge, std::vector<std::size_t> &dof_indices) const = 0;
};

/*! @brief Simple dof map, which orders the dofs by the formula
 *          edge_id * num_components * num_basis_functions + num_basis_functions * component_id + basis_function_id
 */
class SimpleDofMapNetwork : public DofMapNetwork {
public:
  SimpleDofMapNetwork(std::size_t num_components, std::size_t num_basis_functions);

  void dof_indices(const Edge &edge, std::vector<std::size_t> &dof_indices) const override;

private:
  std::size_t d_num_components;
  std::size_t d_num_basis_functions;
};

} // namespace macrocirculation

#endif //TUMORMODELS_DOF_MAP_NETWORK_HPP
