////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tip_vertex_dof_integrator.hpp"

#include <iostream>
#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"

namespace macrocirculation {

TipVertexDofIntegrator::TipVertexDofIntegrator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_total_integration_time(0) {
  reset();
}

void TipVertexDofIntegrator::reset() {
  d_total_integration_time = 0;

  d_quantities.clear();

  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto &v = *d_graph->get_vertex(v_id);
    if (d_graph->get_vertex(v_id)->is_vessel_tree_outflow()) {
      auto &ldofmap = d_dof_map->get_local_dof_map(v);
      d_quantities[v_id] = std::vector<double>(ldofmap.num_local_dof(), 0);
    }
  }
}

double TipVertexDofIntegrator::get_integration_time() const {
  return d_total_integration_time;
}

void TipVertexDofIntegrator::update_vertex_dof(const PetscVec &u, double tau) {
  d_total_integration_time += tau;
  std::vector<double> vertex_dof_values;
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_vessel_tree_outflow()) {
      // get current vertex dof values
      auto local_dof_map = d_dof_map->get_local_dof_map(*v);
      const auto &dof_indices = local_dof_map.dof_indices();
      vertex_dof_values.resize(dof_indices.size());
      extract_dof(dof_indices, u, vertex_dof_values);
      // add values to integral;
      auto &old = d_quantities.at(v_id);
      assert(old.size() == vertex_dof_values.size());
      for (size_t k = 0; k < old.size(); k += 1)
        old[k] += vertex_dof_values[k] * tau;
    }
  }
}

std::map<size_t, std::vector<double>> TipVertexDofIntegrator::get_integral_value(const std::vector<size_t> &vertex_dof_numbers) const {
  std::map<size_t, std::vector<double>> data;
  for (auto v_id : d_graph->get_vertex_ids()) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_vessel_tree_outflow()) {
      data[v_id] = std::vector<double>(vertex_dof_numbers.size(), 0);
      auto &e = *d_graph->get_edge(v->get_edge_neighbors()[0]);
      if (e.rank() == mpi::rank(MPI_COMM_WORLD)) {
        for (size_t k = 0; k < vertex_dof_numbers.size(); k += 1)
          data.at(v_id).at(k) = d_quantities.at(v_id)[k];
      }
      MPI_Bcast(&data[v_id].front(), vertex_dof_numbers.size(), MPI_DOUBLE, e.rank(), d_comm);
    }
  }
  return data;
}

} // namespace macrocirculation