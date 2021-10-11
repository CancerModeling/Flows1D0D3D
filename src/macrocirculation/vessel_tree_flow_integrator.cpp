////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "vessel_tree_flow_integrator.hpp"
#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include <iostream>

namespace macrocirculation {


VesselTreeFlowIntegrator::VesselTreeFlowIntegrator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map)
    : d_comm(comm), d_graph(std::move(graph)), d_dof_map(std::move(dof_map)) {
  reset();
}

void VesselTreeFlowIntegrator::reset() {
  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto vertex = d_graph->get_vertex(v_id);
    if (vertex->is_vessel_tree_outflow()) {
      d_avg_data[v_id] = {0, 0, 0, 0};
    }
  }
}

void VesselTreeFlowIntegrator::add(const PetscVec &u, double tau) {
  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto vertex = d_graph->get_vertex(v_id);
    if (vertex->is_vessel_tree_outflow()) {
      auto &dof_indices = d_dof_map->get_local_dof_map(*vertex).dof_indices();
      std::vector<double> dof_values(dof_indices.size());
      extract_dof(dof_indices, u, dof_values);
      auto p_1d = dof_values.back();
      auto &data = vertex->get_vessel_tree_data();
      auto p_3d = data.p_out;
      auto R = data.resistances.back();
      auto level = data.resistances.size();
      auto flow = std::pow(2, level - 1) / R * (p_1d - p_3d);
      d_avg_data[v_id].flow += tau * flow;
      d_avg_data[v_id].pressure_3d += tau * p_3d * 1e3;
      d_avg_data[v_id].pressure_1d += tau * p_1d * 1e3;
      d_avg_data[v_id].time += tau;
    }
  }
}

std::vector<VesselTreeFlowIntegratorResult> VesselTreeFlowIntegrator::calculate() {
  std::vector<VesselTreeFlowIntegratorResult> results;
  for (const auto &v_id : d_graph->get_vertex_ids()) {
    auto vertex = d_graph->get_vertex(v_id);
    if (!vertex->is_leaf() || !vertex->is_vessel_tree_outflow())
      continue;
    auto edge = d_graph->get_edge(vertex->get_edge_neighbors().front());
    bool owned = edge->rank() == mpi::rank(d_comm);
    std::array<double, 3> data = {0, 0, 0};
    if (owned) {
      double time = d_avg_data.at(v_id).time;
      data[0] = d_avg_data.at(v_id).flow / time;
      data[1] = d_avg_data.at(v_id).pressure_1d / time;
      data[2] = d_avg_data.at(v_id).pressure_3d / time;
    }
    MPI_Bcast(&data[0], data.size(), MPI_DOUBLE, edge->rank(), d_comm);
    if(!edge->has_embedding_data())
      throw std::runtime_error("VesselTreeFlowIntegrator::calculate: needs point coordinates");
    auto p = edge->is_pointing_to(v_id) ? edge->get_embedding_data().points.back() : edge->get_embedding_data().points.front();
    results.push_back({p, v_id, vertex->get_vessel_tree_data().resistances.size(), data[0], data[2]});
  }
  return results;
}

void VesselTreeFlowIntegrator::write(const std::string &folder, const std::string &dataset_name) {
  auto results = calculate();

  if (mpi::rank(d_comm) == 0) {
    write(folder, dataset_name, *d_graph, results);
  }
}

void VesselTreeFlowIntegrator::write(const std::string &folder, const std::string &dataset_name, const GraphStorage &graph, const std::vector<VesselTreeFlowIntegratorResult> &results) {
  using json = nlohmann::json;

  json j;

  auto vertices_list = json::array();

  size_t index = 0;

  for (const auto &v_id : graph.get_vertex_ids()) {
    auto vertex = graph.get_vertex(v_id);
    if (!vertex->is_leaf() || !vertex->is_vessel_tree_outflow())
      continue;

    auto point = results[index].point;

    json vessel_obj = {
      {"vertex_id", v_id},
      {"name", vertex->get_name()},
      {"neighbor_edge_id", vertex->get_edge_neighbors()[0]},
      {"point", {point.x, point.y, point.z}},
      {"average_flow", results[index].averaged_flow},
      // {"average_3d_pressure", results[index].averaged_3D_pressure},
    };

    vertices_list.push_back(vessel_obj);
    index += 1;
  }

  j["vertices"] = vertices_list;

  std::fstream file(folder + "/" + dataset_name + ".json", std::ios::out);
  file << j.dump(1);
}

} // namespace macrocirculation
