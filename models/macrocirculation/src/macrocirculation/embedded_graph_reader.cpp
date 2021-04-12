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
  json j;
  file >> j;

  const std::size_t num_vertices = j["num_vertices"];

  // create the vertices for the given mesh
  std::vector<Vertex *> vertices;
  for (std::size_t idx = 0; idx < num_vertices; idx += 1)
    vertices.push_back(graph.create_vertex().get());

  // create the edges between the vertices
  for (auto vessel : j["vessels"]) {
    size_t left_vertex_id = vessel["left_vertex_id"];
    size_t right_vertex_id = vessel["right_vertex_id"];

    size_t num_micro_edges = vessel["abstract_coordinates"].size() - 1;

    auto edge = graph.connect(*vertices[left_vertex_id], *vertices[right_vertex_id], num_micro_edges);

    std::vector<Point> points;
    for (auto p : vessel["embedded_coordinates"])
      points.emplace_back(p[0], p[1], p[2]);

    edge->add_embedding_data({ points });

    double r_avg = 0;
    for (double r : vessel["radii"])
      r_avg += r;
    r_avg /= vessel["radii"].size();

    const double wall_thickness = vessel["wall_thickness"];
    double E = vessel["elastic_modulus"];
    E /= 100.;

    const double A0 = std::pow(r_avg, 2) * M_PI;
    // const double G0 = calculate_G0(d_wall_width, d_elastic_modulus, d_poisson_ratio, A0);
    const double G0 = 4.0 / 3.0 * std::sqrt(M_PI) * E * wall_thickness / std::sqrt(A0);
    const double length = vessel["vessel_length"];

    edge->add_physical_data({
      G0,
      A0,
      d_rho,
      length
    });
  }
}

void EmbeddedGraphReader::set_parameter(double poisson_ratio, double wall_width, double nu, double rho) {
  d_elastic_modulus = poisson_ratio;
  d_wall_width = wall_width;
  d_poisson_ratio = nu;
  d_rho = rho;
}

} // namespace macrocirculation
