////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "embedded_graph_reader.hpp"
#include "graph_storage.hpp"
#include "vessel_formulas.hpp"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

namespace macrocirculation {

void EmbeddedGraphReader::append(const std::string &filepath, GraphStorage &graph) const {
  std::cout << "WARNING: appending graph with avg(r) and artificial E and nu" << std::endl;

  using json = nlohmann::json;

  std::fstream file(filepath, std::ios::in);
  if (!file.good())
    throw std::runtime_error("file " + filepath + " could not be opened");

  json j;
  file >> j;

  const std::size_t num_vertices = j["vertices"].size();

  const size_t vertex_offset = graph.num_vertices();

  // create the vertices for the given mesh
  for (std::size_t idx = 0; idx < num_vertices; idx += 1)
    graph.create_vertex();

  // create the edges between the vertices
  for (auto vessel : j["vessels"]) {
    size_t left_vertex_id = vertex_offset + vessel["left_vertex_id"].get<size_t>();
    size_t right_vertex_id = vertex_offset + vessel["right_vertex_id"].get<size_t>();

    assert(left_vertex_id < num_vertices + vertex_offset);
    assert(right_vertex_id < num_vertices + vertex_offset);

    size_t num_micro_edges = 0;
    if (vessel.contains("number_edges"))
      num_micro_edges = vessel["number_edges"];
    else if (vessel.contains("abstract_coordinates"))
      num_micro_edges = vessel["abstract_coordinates"].size() - 1;
    else
      throw std::runtime_error("cannot infer number of micro edges");

    auto edge = graph.connect(*graph.get_vertex(left_vertex_id), *graph.get_vertex(right_vertex_id), num_micro_edges);

    if (vessel.contains("embedded_coordinates")) {
      std::vector<Point> points;
      for (auto p : vessel["embedded_coordinates"])
        points.emplace_back(p[0], p[1], p[2]);
      edge->add_embedding_data({points});
    }

    double r_avg = 0;
    if (vessel.contains("radius")) {
      r_avg = vessel["radius"];
    } else if (vessel.contains("radii")) {
      for (double r : vessel["radii"])
        r_avg += r;
      r_avg /= static_cast<double>(vessel["radii"].size());
    } else {
      throw std::runtime_error("cannot infer radius");
    }

    if (vessel.contains("name"))
      edge->set_name(vessel["name"]);

    edge->add_physical_data(PhysicalData::set_from_data(vessel["elastic_modulus"], vessel["wall_thickness"], d_rho, vessel["gamma"], r_avg, vessel["vessel_length"]));
  }

  for (auto vertex : j["vertices"]) {
    size_t id = vertex["id"];

    auto &v = *graph.get_vertex(id + vertex_offset);

    if (!v.is_leaf() && (vertex.contains("peripheral_resistance") || vertex.contains("peripheral_compliance")))
      throw std::runtime_error("malformed network");

    // if the vertex has a name we add it, otherwise we generate it.
    if (vertex.contains("name")) {
      v.set_name(vertex["name"]);
    } else {
      std::stringstream name;
      name << "bifurcation";
      for (auto e_id : v.get_edge_neighbors())
        name << "_" << e_id;
      v.set_name(name.str());
    }

    if (v.is_leaf()) {
      if (vertex.contains("peripheral_resistance") || vertex.contains("peripheral_compliance")) {
        // if (false) {
        const double r = vertex["peripheral_resistance"];
        const double c = vertex["peripheral_compliance"];
        v.set_to_windkessel_outflow(r, c);
      } else {
        v.set_to_free_outflow();
      }
    }
  }
}

void EmbeddedGraphReader::set_boundary_data(const std::string &filepath, GraphStorage &graph) const {
  using json = nlohmann::json;

  std::fstream file(filepath, std::ios::in);
  if (!file.good())
    throw std::runtime_error("file " + filepath + " could not be opened");

  json j;
  file >> j;

  for (const auto &vertex : j["vertices"]) {
    if (!vertex.contains("name"))
      throw std::runtime_error("boundary data can only be set for named vertices");

    auto &v = *graph.find_vertex_by_name(vertex["name"]);

    if (v.is_leaf()) {
      if (vertex.contains("peripheral_resistance") || vertex.contains("peripheral_compliance")) {
        // if (false) {
        const double r = vertex["peripheral_resistance"];
        const double c = vertex["peripheral_compliance"];
        v.set_to_windkessel_outflow(r, c);
      } else {
        v.set_to_free_outflow();
      }
    }
  }
}

std::vector<InputPressuresResults> read_input_pressures(const std::string &filepath) {
  using json = nlohmann::json;

  std::fstream file(filepath, std::ios::in);
  json j;
  file >> j;

  std::vector< InputPressuresResults > results;
  for (auto d : j["vertices"])
  {
    InputPressuresResults result;
    result.name = d["name"].get<std::string>();
    result.t = d["t"].get<std::vector<double>>();
    result.p = d["p"].get<std::vector<double>>();
    result.periodic = d["periodic"].get<bool>();
    results.push_back(result);
  }

  return results;
}

} // namespace macrocirculation
