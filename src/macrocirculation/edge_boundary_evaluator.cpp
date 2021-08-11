////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "edge_boundary_evaluator.hpp"

#include <cmath>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

EdgeBoundaryEvaluator::EdgeBoundaryEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t component)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_component(component),
      d_edge_boundary_communicator(Communicator::create_edge_boundary_value_communicator(comm, d_graph)),
      d_macro_edge_boundary_value(2 * d_graph->num_edges(), NAN) {}

void EdgeBoundaryEvaluator::init(const std::vector<double> &u_prev) {
  evaluate_macro_edge_boundary_values(u_prev);
  d_edge_boundary_communicator.update_ghost_layer(d_macro_edge_boundary_value);
}

void EdgeBoundaryEvaluator::init(const PetscVec &u_prev) {
  evaluate_macro_edge_boundary_values(u_prev);
  d_edge_boundary_communicator.update_ghost_layer(d_macro_edge_boundary_value);
}

template<typename VectorType>
void EdgeBoundaryEvaluator::evaluate_macro_edge_boundary_values(const VectorType &u_prev) {
  std::vector<std::size_t> dof_indices(4, 0);
  std::vector<double> local_dofs(4, 0);

  // as a precaution we fill the boundary value vector with NANs.
  std::fill(d_macro_edge_boundary_value.begin(), d_macro_edge_boundary_value.end(), NAN);

  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto edge = d_graph->get_edge(e_id);
    const auto &param = edge->get_physical_data();
    const auto &local_dof_map = d_dof_map->get_local_dof_map(*edge);
    const double h = param.length / static_cast<double>(local_dof_map.num_micro_edges());

    FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    dof_indices.resize(local_dof_map.num_basis_functions());
    local_dofs.resize(local_dof_map.num_basis_functions());

    local_dof_map.dof_indices(0, d_component, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_macro_edge_boundary_value[2 * edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

    local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, d_component, dof_indices);
    extract_dof(dof_indices, u_prev, local_dofs);
    d_macro_edge_boundary_value[2 * edge->get_id() + 1] = fe.evaluate_dof_at_boundary_points(local_dofs).right;
  }
}

void EdgeBoundaryEvaluator::operator()(const Edge &edge, std::vector<double> &values) const {
  values[0] = d_macro_edge_boundary_value[2 * edge.get_id() + 0];
  values[1] = d_macro_edge_boundary_value[2 * edge.get_id() + 1];
}

double EdgeBoundaryEvaluator::operator()(const Vertex &vertex, const Edge &edge) const {
  if (edge.is_pointing_to(vertex.get_id()))
    return d_macro_edge_boundary_value[2 * edge.get_id() + 1];
  else
    return d_macro_edge_boundary_value[2 * edge.get_id() + 0];
}

void EdgeBoundaryEvaluator::operator()(const Vertex &vertex, std::vector<double> &values) const {
  values.resize(vertex.get_edge_neighbors().size());
  for (size_t k = 0; k < vertex.get_edge_neighbors().size(); k += 1) {
    auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[k]);
    values[k] = (*this)(vertex, edge);
  }
}

} // namespace macrocirculation