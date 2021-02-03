////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_storage.hpp"
#include <algorithm>

namespace macrocirculation {

Primitive::Primitive(std::size_t id) : p_id(id) {}

std::size_t Primitive::get_id() const { return p_id; }

Vertex::Vertex(std::size_t id, const Point &coordinate)
    : Primitive(id),
      p_coordinate(coordinate) {}

const Point &Vertex::get_coordinate() const { return p_coordinate; }

const std::vector<std::size_t> &Vertex::get_edge_neighbors() const { return p_neighbors; }

bool Vertex::is_leaf() const { return p_neighbors.size() == 1; };

bool Vertex::is_unconnected() const { return p_neighbors.empty(); }

Edge::Edge(std::size_t id, std::size_t type_id, const Vertex &v1, const Vertex &v2)
    : Primitive(id),
      p_type_id(type_id),
      p_neighbors({v1.get_id(), v2.get_id()}),
      p_coordinates({v1.get_coordinate(), v2.get_coordinate()}){};

const std::vector<std::size_t> &Edge::get_vertex_neighbors() const { return p_neighbors; };

const Point &Edge::get_coordinate_v0() const { return p_coordinates[0]; };

const Point &Edge::get_coordinate_v1() const { return p_coordinates[1]; };

double Edge::get_length() const {
  return (p_coordinates[0] - p_coordinates[1]).norm();
}

Point Edge::get_midpoint() const {
  return 0.5 * (p_coordinates[0] - p_coordinates[1]);
}

GraphStorage::GraphStorage()
    : p_next_edge_id(0), p_next_vertex_id(0){};

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

std::shared_ptr<Vertex> GraphStorage::create_vertex(const Point &coordinates) {
  const auto id = p_next_vertex_id++;
  auto vertex = std::make_shared<Vertex>(id, coordinates);
  p_vertices[id] = vertex;
  return vertex;
};

std::shared_ptr<Edge> GraphStorage::connect(Vertex &v1, Vertex &v2, std::size_t edge_type_id) {
  return connect(v1, v2, edge_type_id, p_next_edge_id++);
}

std::shared_ptr<Edge> GraphStorage::connect(Vertex &v1, Vertex &v2, std::size_t edge_type_id, std::size_t edge_id) {
  if (p_vertices.find(v1.get_id()) == p_vertices.end() || p_vertices.find(v2.get_id()) == p_vertices.end())
    throw std::runtime_error("vertices not found in storage");

  if (v1.get_id() == v2.get_id())
    throw std::runtime_error("connecting same vertex");

  if (p_edges.find(edge_id) != p_edges.end())
    throw std::runtime_error("edge with given id already in the graph");

  auto edge = std::make_shared<Edge>(edge_id, edge_type_id, v1, v2);
  p_edges[edge_id] = edge;

  v1.p_neighbors.push_back(edge->get_id());
  v2.p_neighbors.push_back(edge->get_id());

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

RefinementData GraphStorage::refine(Edge &e) {
  auto v0 = get_vertex(e.get_vertex_neighbors()[0]);
  auto v1 = get_vertex(e.get_vertex_neighbors()[1]);
  auto v0p5 = create_vertex(0.5 * (v0->get_coordinate() + v1->get_coordinate()));

  auto id = e.get_id();

  remove(e);

  auto e0 = connect(*v0, *v0p5, e.p_type_id, id);
  auto e1 = connect(*v0p5, *v1, e.p_type_id);

  return RefinementData(v0p5, e0, e1);
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

void GraphStorage::refine(std::size_t num_refinements) {
  for (std::size_t ridx = 0; ridx < num_refinements; ridx += 1) {
    for (auto &id : get_edge_ids()) {
      auto edge = get_edge(id);
      refine(*edge);
    }
  }
}

double GraphStorage::num_edges() const {
  return p_edges.size();
}

void GraphDataWriter::add_midpoint_data(const std::string &name, std::vector<double> data)
{
  midpoint_data.push_back(std::move(NamedField(name, std::move(data))));
}

std::string GraphDataWriter::get_path(const std::string & filename, std::size_t time_idx) const
{
  std::string path = filename + "_";
  path += std::to_string(time_idx);
  path += ".vtk";
  return path;
}

void GraphDataWriter::write_vtk(const std::string &filename, GraphStorage &storage, std::size_t time_idx) const
{
  const auto path = get_path(filename, time_idx);

  const std::size_t number_of_points = storage.num_edges() * 2;
  const std::size_t number_of_segments = storage.num_edges();

  std::fstream filevtk;
  filevtk.open(path, std::ios::out);
  filevtk << "# vtk DataFile Version 2.0" << std::endl;
  filevtk << "Network Nutrient Transport" << std::endl;
  filevtk << "ASCII" << std::endl;
  filevtk << "DATASET POLYDATA" << std::endl;

  filevtk << "POINTS " << number_of_points << " float" << std::endl;
  // the part for our point data
  for (auto e_id : storage.get_edge_ids()) {
    auto edge = storage.get_edge(e_id);
    const auto &p0 = edge->get_coordinate_v0();
    const auto &p1 = edge->get_coordinate_v1();

    // write the two points into our file
    filevtk << p0(0) << " " << p0(1) << " " << p0(2) << std::endl;
    filevtk << p1(0) << " " << p1(1) << " " << p1(2) << std::endl;
  }
  filevtk << " " << std::endl;

  filevtk << "LINES " << number_of_segments << " " << 3 * number_of_segments << std::endl;
  for (std::size_t idx=0; idx< number_of_segments; idx+=1)
    filevtk << "2 " << 2*idx << " " << 2*idx+1 << std::endl;
  filevtk << " " << std::endl;

  filevtk << "CELL_DATA " << number_of_segments << std::endl;
  for (std::size_t idx = 0; idx < midpoint_data.size(); idx += 1) {
    filevtk << "SCALARS " << midpoint_data[idx].name << " float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    for (auto e_id : storage.get_edge_ids()) {
      auto edge = storage.get_edge(e_id);
      filevtk << midpoint_data[idx].values[edge->get_id()] << std::endl;
    }
  }
}


} // namespace macrocirculation
