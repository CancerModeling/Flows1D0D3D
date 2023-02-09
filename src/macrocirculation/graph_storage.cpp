////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#include "graph_storage.hpp"
#include <cmath> // (add before graph_storage.hpp) easy fix for compile error about std::abs

#include <algorithm>
#include <cassert>
#include <utility>

#include "vessel_formulas.hpp"

namespace macrocirculation {

Point::Point(double x, double y, double z)
    : x(x), y(y), z(z) {}

Point convex_combination(const Point &left, const Point &right, double theta) {
  return {(1 - theta) * left.x + theta * right.x,
          (1 - theta) * left.y + theta * right.y,
          (1 - theta) * left.z + theta * right.z};
}

PhysicalData::PhysicalData(double elastic_modulus, double G0, double A0, double rho, double length, double viscosity, double gamma, double radius)
    : elastic_modulus(elastic_modulus),
      G0(G0),
      A0(A0),
      rho(rho),
      length(length),
      viscosity(viscosity),
      gamma(gamma),
      radius(radius) {}

std::size_t DiscretizationData::num_micro_edges() const {
  return lengths.size();
}

Primitive::Primitive(std::size_t id)
    : p_id(id) {}

std::size_t Primitive::get_id() const {
  return p_id;
}

PhysicalData PhysicalData::set_from_data(double elastic_modulus, double wall_thickness, double density, double gamma, double radius, double length) {
  // [E] = Pa = N m^{-2} = kg s^{-2} m^{-1} = kg s^{-2} cm^{-1} 100^{-1}
  const double E = elastic_modulus / 100;

  const double A0 = std::pow(radius, 2) * M_PI;

  const double G0 = calculate_G0(wall_thickness, E, A0);

  // the viscosity
  const double viscosity = viscosity_bloodplasma(radius);

  return {elastic_modulus, G0, A0, density, length, viscosity, gamma, radius};
}

double PhysicalData::get_c0() const {
  return std::pow(G0 / (2.0 * rho), 0.5);
}

InterGraphConnection::InterGraphConnection(const std::shared_ptr<GraphStorage> &graph, const Vertex &vertex)
    : d_graph(graph),
      d_vertex_id(vertex.get_id()) {}

const GraphStorage &InterGraphConnection::get_graph() const {
  auto graph = d_graph.lock();
  return *graph;
}

size_t Vertex::local_edge_index(const Edge &edge) const {
  for (size_t k = 0; k < get_edge_neighbors().size(); k += 1)
    if (get_edge_neighbors()[k] == edge.get_id())
      return k;
  throw std::runtime_error("vertex (vid=" + std::to_string(get_id()) + ") does not have given edge neighbor (eid=" + std::to_string(edge.get_id()));
}

const Vertex &InterGraphConnection::get_vertex() const {
  return *get_graph().get_vertex(d_vertex_id);
}


const std::string &Primitive::get_name() const {
  return p_name;
}

void Primitive::set_name(const std::string &name) {
  p_name = name;
}

double default_inflow_function(double) {
  throw std::runtime_error("inflow value at inflow boundary not set");
}

Vertex::Vertex(std::size_t id)
    : Primitive(id),
      p_inflow_value(default_inflow_function),
      p_flow_type(FlowType::Undefined),
      d_bcs_finalized(false) {}

const std::vector<std::size_t> &Vertex::get_edge_neighbors() const {
  return p_neighbors;
}

bool Vertex::is_leaf() const {
  return p_neighbors.size() == 1;
};

bool Vertex::is_unconnected() const {
  return p_neighbors.empty();
}

bool Vertex::is_bifurcation() const {
  // TODO: move quick and dirty sanity check somewhere else
  assert(is_leaf() ? p_flow_type != FlowType::Undefined : true);
  return p_neighbors.size() > 1;
}

void Vertex::set_to_inflow_with_fixed_flow(std::function<double(double)> inflow_value) {
  if (!is_leaf())
    throw std::runtime_error("inflow bc can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::InflowFixedFlow;
  p_inflow_value = std::move(inflow_value);
}

void Vertex::set_to_inflow_with_fixed_pressure(std::function<double(double)> pressure_value) {
  if (!is_leaf())
    throw std::runtime_error("inflow bc can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::InflowFixedPressure;
  p_inflow_value = std::move(pressure_value);
}

void Vertex::set_to_free_outflow() {
  if (!is_leaf())
    throw std::runtime_error("free outflow bc can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::FreeOutflow;
  p_inflow_value = default_inflow_function;
}

void Vertex::set_to_windkessel_outflow(double r, double c) {
  if (!is_leaf())
    throw std::runtime_error("windkessel bc can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::Windkessel;
  p_peripheral_vessel_data.resistance = r;
  p_peripheral_vessel_data.compliance = c;
  // TODO: Make this settable from outside!
  p_peripheral_vessel_data.p_out = 5.0 * 1.333322;
}

void Vertex::set_to_vessel_tree_outflow(double p, const std::vector<double> &resistances, const std::vector<double> &capacitances, const std::vector<double> &radii, size_t furcation_number) {
  if (!is_leaf())
    throw std::runtime_error("tree bc can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::VesselTree;
  p_vessel_tree_data.p_out = p;
  p_vessel_tree_data.resistances = resistances;
  p_vessel_tree_data.capacitances = capacitances;
  p_vessel_tree_data.radii = radii;
  p_vessel_tree_data.furcation_number = furcation_number;
}

bool Vertex::bc_finalized() const { return d_bcs_finalized; }

void Vertex::finalize_bcs() { d_bcs_finalized = true; }

void Vertex::update_vessel_tip_pressures(double p) {
  if (is_vessel_tree_outflow())
    p_vessel_tree_data.p_out = p;
  else if (is_rcl_outflow())
    p_rcl_data.p_out = p;
  else if (is_windkessel_outflow())
    p_peripheral_vessel_data.p_out = p;
  else
    throw std::runtime_error("updating the boundary pressure was not implemented for the given boundary type.");
}

const PeripheralVesselData &Vertex::get_peripheral_vessel_data() const {
  assert(is_windkessel_outflow());
  return p_peripheral_vessel_data;
}

const VesselTreeData &Vertex::get_vessel_tree_data() const {
  assert(is_vessel_tree_outflow());
  return p_vessel_tree_data;
}

const RCLModel &Vertex::get_rcl_data() const {
  assert(is_rcl_outflow());
  return p_rcl_data;
}

bool Vertex::is_rcl_outflow() const {
  return p_flow_type == FlowType::RCLModel;
}

const LinearCharacteristicData &Vertex::get_linear_characteristic_data() const {
  assert(is_linear_characteristic_inflow());
  return p_linear_characteristic_data;
}

const NonlinearCharacteristicData &Vertex::get_nonlinear_characteristic_data() const {
  assert(is_nonlinear_characteristic_inflow());
  return p_nonlinear_characteristic_data;
}

bool Vertex::is_free_outflow() const { return p_flow_type == FlowType::FreeOutflow; }

bool Vertex::is_windkessel_outflow() const { return p_flow_type == FlowType::Windkessel; };

bool Vertex::is_vessel_tree_outflow() const { return p_flow_type == FlowType::VesselTree; };

bool Vertex::is_inflow_with_fixed_flow() const {
  return p_flow_type == FlowType::InflowFixedFlow;
}

bool Vertex::is_inflow_with_fixed_pressure() const {
  return p_flow_type == FlowType::InflowFixedPressure;
}

bool Vertex::is_linear_characteristic_inflow() const {
  return p_flow_type == FlowType::LinearCharacteristic;
}

bool Vertex::is_nonlinear_characteristic_inflow() const {
  return p_flow_type == FlowType::NonlinearCharacteristic;
}

double Vertex::get_inflow_value(double time) const {
  return p_inflow_value(time);
}

void Vertex::set_to_linear_characteristic_inflow(double C, double L, bool points_towards_vertex, double p, double q) {
  if (!is_leaf())
    throw std::runtime_error("linear characteristic inflow can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::LinearCharacteristic;
  double sigma = points_towards_vertex ? +1 : -1;
  p_linear_characteristic_data.C = C;
  p_linear_characteristic_data.L = L;
  p_linear_characteristic_data.points_towards_vertex = points_towards_vertex;
  p_linear_characteristic_data.p = p;
  p_linear_characteristic_data.q = sigma * q;
}

void Vertex::set_to_vessel_rcl_outflow(double p, const std::vector<double> &resistances, const std::vector<double> &capacitances, const std::vector<double> &inductances) {
  if (!is_leaf())
    throw std::runtime_error("rcl outflow can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::RCLModel;
  p_rcl_data.resistances = resistances;
  p_rcl_data.capacitances = capacitances;
  p_rcl_data.inductances = inductances;
  p_rcl_data.p_out = p;
}

void Vertex::set_to_nonlinear_characteristic_inflow(double G0, double A0, double rho, bool points_towards_vertex, double p, double q) {
  if (!is_leaf())
    throw std::runtime_error("nonlinear characteristic inflow can only be set for leaf nodes (vertex name = " + get_name() + ")");
  if (d_bcs_finalized)
    throw std::runtime_error("finalized boundary conditions cannot be changed.");
  p_flow_type = FlowType::NonlinearCharacteristic;
  double sigma = points_towards_vertex ? +1 : -1;
  p_nonlinear_characteristic_data.G0 = G0;
  p_nonlinear_characteristic_data.A0 = A0;
  p_nonlinear_characteristic_data.rho = rho;
  p_nonlinear_characteristic_data.points_towards_vertex = points_towards_vertex;
  p_nonlinear_characteristic_data.p = p;
  p_nonlinear_characteristic_data.q = sigma * q;
}

void Vertex::update_linear_characteristic_inflow(double p, double q) {
  assert(p_flow_type == FlowType::LinearCharacteristic);
  double sigma = p_linear_characteristic_data.points_towards_vertex ? +1 : -1;
  p_linear_characteristic_data.p = p;
  p_linear_characteristic_data.q = sigma * q;
}

void Vertex::update_nonlinear_characteristic_inflow(double p, double q) {
  assert(p_flow_type == FlowType::NonlinearCharacteristic);
  double sigma = p_nonlinear_characteristic_data.points_towards_vertex ? +1 : -1;
  p_nonlinear_characteristic_data.p = p;
  p_nonlinear_characteristic_data.q = sigma * q;
}

void Vertex::add_inter_graph_connection(std::shared_ptr<GraphStorage> graph, Vertex &v) {
  d_inter_graph_connections.emplace_back(graph, v);
}

const std::vector<InterGraphConnection> &Vertex::get_inter_graph_connections() const {
  return d_inter_graph_connections;
}

void Vertex::connect(const std::shared_ptr<GraphStorage> &g1, Vertex &v1, const std::shared_ptr<GraphStorage> &g2, Vertex &v2) {
  v1.add_inter_graph_connection(g2, v2);
  v2.add_inter_graph_connection(g1, v1);
}

bool Edge::is_pointing_to(std::size_t vertex_id) const {
  return p_neighbors[1] == vertex_id;
}

const PhysicalData &Edge::get_physical_data() const {
  return *physical_data;
}

const DiscretizationData &Edge::get_discretization_data() const {
  return *discretization_data;
}

const EmbeddingData &Edge::get_embedding_data() const {
  return *embedding_data;
}

void Edge::add_physical_data(const PhysicalData &data) {
  physical_data = std::make_unique<PhysicalData>(data);
}

void Edge::add_discretization_data(const DiscretizationData &data) {
  discretization_data = std::make_unique<DiscretizationData>(data);
};

void Edge::add_embedding_data(const EmbeddingData &data) {
  embedding_data = std::make_unique<EmbeddingData>(data);
};

bool Edge::has_physical_data() const {
  return physical_data != nullptr;
}

bool Edge::has_discretization_data() const {
  return discretization_data != nullptr;
}

bool Edge::has_embedding_data() const {
  return embedding_data != nullptr;
}

void Edge::assign_to_rank(int rank) { d_rank = rank; }

int Edge::rank() const { return d_rank; }

void GraphStorage::reorder_edges(Vertex &v) {
  if (v.get_edge_neighbors().size() > 2)
    throw std::runtime_error("ordering edges on vertex with more than two neighbors.");

  // nothing to reorder
  if (v.get_edge_neighbors().size() <= 1)
    return;

  const auto e0 = get_edge(v.get_edge_neighbors()[0]);
  const auto e1 = get_edge(v.get_edge_neighbors()[1]);

  if (e0->is_pointing_to(v.get_id()) == e1->is_pointing_to(v.get_id()))
    throw std::runtime_error("edges cannot be reordered");

  // nothing to reorder
  if (e0->is_pointing_to(v.get_id()))
    return;

  std::swap(v.p_neighbors[0], v.p_neighbors[1]);
}

Edge::Edge(std::size_t id,
           const Vertex &v1,
           const Vertex &v2,
           std::size_t first_micro_edge_id,
           std::size_t first_micro_vertex_id,
           std::size_t num_micro_edges)
    : Primitive(id),
      p_neighbors({v1.get_id(), v2.get_id()}),
      d_rank(0),
      d_micro_edges(),
      d_micro_vertices() {
  assert(num_micro_edges > 0);

  // create micro edges
  for (std::size_t local_micro_edge_id = 0; local_micro_edge_id < num_micro_edges; local_micro_edge_id += 1)
    d_micro_edges.emplace_back(local_micro_edge_id, first_micro_edge_id + local_micro_edge_id);

  // create micro vertices
  for (std::size_t local_micro_vertex_id = 0; local_micro_vertex_id < num_micro_edges + 1; local_micro_vertex_id += 1)
    d_micro_vertices.emplace_back(local_micro_vertex_id, first_micro_vertex_id + local_micro_vertex_id);

  // connect inner micro_vertices
  d_micro_vertices[0].d_right_edge = &d_micro_edges[0];
  d_micro_vertices[0].d_right_vertex = &d_micro_vertices[1];
  for (std::size_t local_micro_vertex_id = 1; local_micro_vertex_id < num_micro_edges; local_micro_vertex_id += 1) {
    d_micro_vertices[local_micro_vertex_id].d_left_edge = &d_micro_edges[local_micro_vertex_id - 1];
    d_micro_vertices[local_micro_vertex_id].d_right_edge = &d_micro_edges[local_micro_vertex_id];
    d_micro_vertices[local_micro_vertex_id].d_left_vertex = &d_micro_vertices[local_micro_vertex_id - 1];
    d_micro_vertices[local_micro_vertex_id].d_right_vertex = &d_micro_vertices[local_micro_vertex_id + 1];
  }
  d_micro_vertices.back().d_left_edge = &d_micro_edges.back();
  d_micro_vertices.back().d_left_vertex = &d_micro_vertices.back() - 1;
};

std::size_t Edge::num_micro_edges() const { return d_micro_edges.size(); };

std::size_t Edge::num_micro_vertices() const { return d_micro_edges.size() + 1; };

const std::vector<MicroEdge> &Edge::micro_edges() const { return d_micro_edges; };

const std::vector<MicroVertex> &Edge::micro_vertices() const { return d_micro_vertices; }

const std::vector<std::size_t> &Edge::get_vertex_neighbors() const {
  return p_neighbors;
};

Edge::InnerVerticesIterator Edge::inner_micro_vertices() const {
  return InnerVerticesIterator{d_micro_vertices};
}

const MicroVertex &Edge::left_micro_vertex() const {
  assert(!d_micro_vertices.empty());
  return d_micro_vertices[0];
}
const MicroVertex &Edge::right_micro_vertex() const {
  assert(!d_micro_vertices.empty());
  return d_micro_vertices.back();
}

std::size_t Edge::get_adajcent_micro_edge_id(const Vertex &vertex) const {
  if (vertex.get_id() != get_vertex_neighbors()[0] && vertex.get_id() != get_vertex_neighbors()[1])
    throw std::runtime_error("Edge::get_adjacent_micro_edge_id: edge is not adjacent to given vertex");

  return is_pointing_to(vertex.get_id()) ? num_micro_edges() - 1 : 0;
}

GraphStorage::GraphStorage()
    : p_next_edge_id(0),
      p_next_vertex_id(0),
      d_num_micro_edges(0),
      d_num_micro_vertices(0){};

std::shared_ptr<Edge> GraphStorage::get_edge(std::size_t id) {
  return p_edges.at(id);
};

std::shared_ptr<const Edge> GraphStorage::get_edge(std::size_t id) const {
  return p_edges.at(id);
}

std::shared_ptr<Vertex> GraphStorage::get_vertex(std::size_t id) {
  return p_vertices.at(id);
}
std::shared_ptr<const Vertex> GraphStorage::get_vertex(std::size_t id) const {
  return p_vertices.at(id);
}

std::shared_ptr<Vertex> GraphStorage::create_vertex() {
  const auto id = p_next_vertex_id++;
  auto vertex = std::make_shared<Vertex>(id);
  p_vertices[id] = vertex;
  return vertex;
};

std::shared_ptr<Edge> GraphStorage::connect(Vertex &v1, Vertex &v2, std::size_t num_micro_edges) {
  return connect(v1, v2, p_next_edge_id++, num_micro_edges);
}

std::shared_ptr<Edge> GraphStorage::connect(Vertex &v1, Vertex &v2, std::size_t edge_id, std::size_t num_local_micro_edges) {
  if (p_vertices.find(v1.get_id()) == p_vertices.end() || p_vertices.find(v2.get_id()) == p_vertices.end())
    throw std::runtime_error("vertices not found in storage");

  if (v1.get_id() == v2.get_id())
    throw std::runtime_error("connecting same vertex");

  if (p_edges.find(edge_id) != p_edges.end())
    throw std::runtime_error("edge with given id already in the graph");

  auto edge = std::make_shared<Edge>(edge_id, v1, v2, d_num_micro_edges, d_num_micro_vertices, num_local_micro_edges);
  p_edges[edge_id] = edge;

  v1.p_neighbors.push_back(edge->get_id());
  v2.p_neighbors.push_back(edge->get_id());

  d_num_micro_edges += num_local_micro_edges;
  d_num_micro_vertices += num_local_micro_edges + 1;

  return edge;
}

void GraphStorage::remove(Edge &e) {
  p_edges.erase(p_edges.find(e.get_id()));
  auto v0 = get_vertex(e.get_vertex_neighbors()[0]);
  auto v1 = get_vertex(e.get_vertex_neighbors()[1]);
  v0->p_neighbors.erase(std::remove(v0->p_neighbors.begin(), v0->p_neighbors.end(), e.get_id()), v0->p_neighbors.end());
  v1->p_neighbors.erase(std::remove(v1->p_neighbors.begin(), v1->p_neighbors.end(), e.get_id()), v1->p_neighbors.end());
  e.p_neighbors.clear();
}

std::vector<std::size_t> GraphStorage::get_edge_ids() const {
  std::vector<std::size_t> keys;
  transform(p_edges.begin(), p_edges.end(), back_inserter(keys), [](auto v) { return v.first; });
  return keys;
}

std::vector<std::size_t> GraphStorage::get_vertex_ids() const {
  std::vector<std::size_t> keys;
  transform(p_vertices.begin(), p_vertices.end(), back_inserter(keys), [](auto v) { return v.first; });
  return keys;
}

size_t GraphStorage::num_vertices() const {
  return p_vertices.size();
}

size_t GraphStorage::num_edges() const {
  return p_edges.size();
}

std::vector<std::size_t> GraphStorage::get_active_edge_ids(int rank) const {
  std::vector<std::size_t> active_edge_ids;
  for (const auto &it : p_edges) {
    auto edge = it.second;

    if (edge->rank() == rank)
      active_edge_ids.push_back(edge->get_id());
  }
  return active_edge_ids;
}

std::vector<std::size_t> GraphStorage::get_active_vertex_ids(int rank) const {
  std::vector<std::size_t> active_vertex_ids;
  for (const auto &v_it : p_vertices) {
    auto vertex = v_it.second;

    if (vertex_is_neighbor_of_rank(*vertex, rank))
      active_vertex_ids.push_back(vertex->get_id());
  }
  return active_vertex_ids;
}

std::vector<std::size_t> GraphStorage::get_active_and_connected_vertex_ids(int rank) const {
  std::vector<std::size_t> active_vertex_ids;
  for (const auto &v_it : p_vertices) {
    auto vertex = v_it.second;

    if (vertex_is_neighbor_of_rank(*vertex, rank) || vertex_is_connected_to_rank(*vertex, rank))
      active_vertex_ids.push_back(vertex->get_id());
  }
  return active_vertex_ids;
}

bool GraphStorage::edge_is_neighbor_of_rank(const Edge &e, int rank) const {
  for (auto v_id : e.get_vertex_neighbors()) {
    auto vertex = get_vertex(v_id);

    if (vertex_is_neighbor_of_rank(*vertex, rank))
      return true;
  }
  return false;
}

bool GraphStorage::vertex_is_connected_to_rank(const Vertex &vertex, int rank) const {
  for (auto &con : vertex.get_inter_graph_connections()) {
    auto &connected_graph = con.get_graph();
    auto &connected_vertex = con.get_vertex();

    if (connected_graph.vertex_is_neighbor_of_rank(connected_vertex, rank))
      return true;
  }
  return false;
}

bool GraphStorage::edge_is_connected_to_rank(const Edge &e, int rank) const {
  for (auto v_id : e.get_vertex_neighbors()) {
    auto vertex = get_vertex(v_id);

    if (vertex_is_connected_to_rank(*vertex, rank))
      return true;
  }
  return false;
}

bool GraphStorage::vertex_is_neighbor_of_rank(const Vertex &v, int rank) const {
  for (const auto &e_id : v.get_edge_neighbors()) {
    auto neighbor_edge = get_edge(e_id);
    if (neighbor_edge->rank() == rank)
      return true;
  }
  return false;
}

std::vector<std::size_t> GraphStorage::get_ghost_edge_ids(int main_rank, int ghost_rank) const {
  std::vector<std::size_t> ghost_edge_ids;
  for (const auto &it : p_edges) {
    auto ghost_rank_edge = it.second;

    // if the edge does not belong to the ghost rank we skip it
    if (ghost_rank_edge->rank() != ghost_rank)
      continue;

    if (edge_is_neighbor_of_rank(*ghost_rank_edge, main_rank) || edge_is_connected_to_rank(*ghost_rank_edge, main_rank))
      ghost_edge_ids.push_back(ghost_rank_edge->get_id());
  }
  return ghost_edge_ids;
}

void GraphStorage::assign_edge_to_rank(Edge &edge, int rank) {
  edge.assign_to_rank(rank);
}

std::vector<std::shared_ptr<Vertex>> GraphStorage::find_embedded_vertices(const Point &p) const {
  std::vector<std::shared_ptr<Vertex>> found_vertices;

  for (const auto &it : p_vertices) {
    const auto vertex = it.second;

    bool has_coordinates = false;
    for (const auto e_id : vertex->get_edge_neighbors()) {
      const auto edge = get_edge(e_id);
      if (edge->has_embedding_data() && edge->get_embedding_data().points.size() >= 2) {
        const auto &points = edge->get_embedding_data().points;
        const auto p_edge = edge->is_pointing_to(vertex->get_id()) ? points.back() : points.front();
        const double norm = std::abs(p_edge.x - p.x) + std::abs(p_edge.y - p.y) + std::abs(p_edge.z - p.z);
        if (norm < 1e-12)
          has_coordinates = true;
      }
    }

    if (has_coordinates)
      found_vertices.push_back(vertex);
  }

  return found_vertices;
}

std::shared_ptr<Edge> GraphStorage::find_edge_by_name(const std::string &name) {
  auto pos = std::find_if(p_edges.begin(), p_edges.end(), [&](auto e) { return e.second->get_name() == name; });

  if (pos == p_edges.end())
    throw std::runtime_error("edge " + name + " not found in graph storage.");

  return pos->second;
}

bool GraphStorage::has_named_vertex(const std::string &name) const {
  auto pos = std::find_if(p_vertices.begin(), p_vertices.end(), [&](auto v) { return v.second->get_name() == name; });
  return pos != p_vertices.end();
}

std::shared_ptr<Vertex> GraphStorage::find_vertex_by_name(const std::string &name) {
  auto pos = std::find_if(p_vertices.begin(), p_vertices.end(), [&](auto v) { return v.second->get_name() == name; });

  if (pos == p_vertices.end())
    throw std::runtime_error("vertex " + name + " not found in graph storage.");

  return pos->second;
}

void GraphStorage::finalize_bcs() {
  for (auto it : p_vertices) {
    auto &vertex = *it.second;
    vertex.finalize_bcs();
  }
}

bool GraphStorage::owns_primitive(const Vertex &vertex, size_t rank) const {
  return vertex_is_neighbor_of_rank(vertex, rank);
}

std::vector<double> get_normals(const GraphStorage &graph, const Vertex &v) {
  std::vector<double> sigma(v.get_edge_neighbors().size(), 0);
  get_normals(graph, v, sigma);
  return sigma;
}

void get_normals(const GraphStorage &graph, const Vertex &v, std::vector<double> &sigma) {
  for (size_t k = 0; k < v.get_edge_neighbors().size(); k += 1) {
    auto &e = *graph.get_edge(v.get_edge_neighbors()[k]);
    sigma[k] = e.is_pointing_to(v.get_id()) ? +1. : -1.;
  }
}

} // namespace macrocirculation
