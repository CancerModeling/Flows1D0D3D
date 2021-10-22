////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "create_3_vessel_network.hpp"

#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

namespace test_macrocirculation {

namespace util {

std::shared_ptr<mc::GraphStorage> create_3_vessel_network() {
  auto graph = std::make_shared<mc::GraphStorage>();
  auto &v0 = *graph->create_vertex();
  auto &v1 = *graph->create_vertex();
  auto &v2 = *graph->create_vertex();
  auto &v3 = *graph->create_vertex();

  v0.set_name("in");
  v1.set_name("out_1");
  v2.set_name("out_2");
  v3.set_name("out_3");

  const double density = 1.028e-3;

  double vessel_length_0 = 4.0;
  double radius_0 = 1.2;
  double wall_thickness_0 = 0.163;
  double elastic_modulus_0 = 400000.0;
  double gamma_0 = 9;
  size_t number_edges_0 = 10;

  double vessel_length_1 = 2.0;
  double radius_1 = 1.12;
  double wall_thickness_1 = 0.126;
  double elastic_modulus_1 = 400000.0;
  double gamma_1 = 9;
  size_t number_edges_1 = 6;

  double vessel_length_2 = 3.4;
  double radius_2 = 0.62;
  double wall_thickness_2 = 0.08;
  double elastic_modulus_2 = 400000.0;
  double gamma_2 = 9;
  size_t number_edges_2 = 10;

  auto data_0 = mc::PhysicalData::set_from_data(elastic_modulus_0, wall_thickness_0, density, gamma_0, radius_0, vessel_length_0);
  auto data_1 = mc::PhysicalData::set_from_data(elastic_modulus_1, wall_thickness_1, density, gamma_1, radius_1, vessel_length_1);
  auto data_2 = mc::PhysicalData::set_from_data(elastic_modulus_2, wall_thickness_2, density, gamma_2, radius_2, vessel_length_2);

  auto &e0 = *graph->connect(v0, v1, number_edges_0);
  e0.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 0, 0), mc::Point(0, 1, 0)}}));
  auto &e1 = *graph->connect(v1, v2, number_edges_1);
  e1.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 1, 0), mc::Point(1, 1, 0)}}));
  auto &e2 = *graph->connect(v1, v3, number_edges_2);
  e2.add_embedding_data(mc::EmbeddingData({{mc::Point(0, 1, 0), mc::Point(-1, 1, 0)}}));

  e0.add_physical_data(data_0);
  e1.add_physical_data(data_1);
  e2.add_physical_data(data_2);

  v0.set_to_inflow_with_fixed_flow(mc::heart_beat_inflow(485.));

  v2.set_to_windkessel_outflow(1.718414143839568, 0.7369003586207183);
  v3.set_to_windkessel_outflow(6.085065767841017, 0.20809964052427588);

  return graph;
}

} // namespace util

} // namespace test_macrocirculation