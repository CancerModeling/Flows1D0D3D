////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "write_0d_data_to_json.hpp"

#include <iomanip>
#include <nlohmann/json.hpp>

#include "communication/mpi.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

void write_0d_data_to_json(const std::string &output_path,
                           const std::shared_ptr<GraphStorage> &graph,
                           const std::vector<double> &t,
                           const std::map<size_t, std::vector<Values0DModel>> &input) {

  // only root writes the file
  if (mpi::rank(MPI_COMM_WORLD) != 0)
    return;

  using json = nlohmann::json;

  json j;

  j["time"] = t;

  j["vertices"] = json::array();

  for (auto it : input) {
    auto v = json::object();

    const auto &vertex = *graph->get_vertex(it.first);

    // we add some information to which vertex this information belongs
    v["name"] = vertex.get_name();
    v["neighbor_edge_id"] = vertex.get_edge_neighbors()[0];
    v["neighbor_edge_name"] = graph->get_edge(vertex.get_edge_neighbors()[0])->get_name();

    std::vector<double> p_c;
    for (auto r : it.second)
      p_c.push_back(r.p_c);
    v["p_c"] = p_c;

    std::vector<double> q;
    for (auto r : it.second)
      q.push_back(r.q);
    v["q"] = q;

    j["vertices"].push_back(v);
  }

  // write back
  std::ofstream o(output_path);
  o << std::setw(1);
  o << j;
}

} // namespace macrocirculation
