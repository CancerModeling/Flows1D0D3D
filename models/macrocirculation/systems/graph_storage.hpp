////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_STORAGE_H
#define TUMORMODELS_GRAPH_STORAGE_H

#include <libmesh/point.h>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <fstream>

namespace macrocirculation {

using libMesh::Point;

class GraphStorage;

class Primitive {
public:
  explicit Primitive(std::size_t id);

  std::size_t get_id() const;

protected:
  std::size_t p_id;
};

class Vertex : public Primitive {
public:
  Vertex(std::size_t id, const Point &coordinates);

  const Point &get_coordinate() const;

  const std::vector<std::size_t> &get_edge_neighbors() const;

  bool is_leaf() const;
  bool is_unconnected() const;

private:
  Point p_coordinate;

  std::vector<std::size_t> p_neighbors;

  friend GraphStorage;
};

class Edge : public Primitive {
public:
  Edge(std::size_t id, std::size_t type_id, const Vertex &v1, const Vertex &v2);

  const std::vector<std::size_t> &get_vertex_neighbors() const;

  const Point &get_coordinate_v0() const;
  const Point &get_coordinate_v1() const;

  double get_length() const;
  Point get_midpoint() const;

protected:
  std::vector<std::size_t> p_neighbors;

  std::vector<Point> p_coordinates;

  std::size_t p_type_id;

  friend GraphStorage;
};

struct RefinementData {
  RefinementData(std::shared_ptr<Vertex> v, std::shared_ptr<Edge> e0, std::shared_ptr<Edge> e1)
      : new_vertex(std::move(v)),
        new_edge_0(std::move(e0)),
        new_edge_1(std::move(e1)) {}

  std::shared_ptr<Vertex> new_vertex;
  std::shared_ptr<Edge> new_edge_0;
  std::shared_ptr<Edge> new_edge_1;
};

class GraphStorage {
public:
  GraphStorage();

  std::shared_ptr<Edge> get_edge(std::size_t id);
  std::shared_ptr<const Edge> get_edge(std::size_t id) const;

  std::shared_ptr<Vertex> get_vertex(std::size_t id);
  std::shared_ptr<const Vertex> get_vertex(std::size_t id) const;

  std::shared_ptr<Vertex> create_vertex(const Point &coordinates);
  std::shared_ptr<Edge> connect(Vertex &v1, Vertex &v2, std::size_t edge_type_id);
  void remove(Edge &e);

  RefinementData refine(Edge &e);
  void refine(std::size_t num_refinements = 1);

  std::vector<std::size_t> get_edge_ids() const;
  std::vector<std::size_t> get_vertex_ids() const;

  double num_edges() const;

private:
  std::shared_ptr<Edge> connect(Vertex &v1, Vertex &v2, std::size_t edge_type_id, std::size_t edge_id);

  std::map<std::size_t, std::shared_ptr<Edge>> p_edges;
  std::map<std::size_t, std::shared_ptr<Vertex>> p_vertices;

  std::size_t p_next_edge_id;
  std::size_t p_next_vertex_id;
};

class GraphDataWriter {
public:
  void add_midpoint_data(const std::string &name, std::vector<double> data);

  void write_vtk(const std::string &filename, GraphStorage &storage, std::size_t time_idx) const;

private:
  struct NamedField {
    NamedField(std::string n, std::vector<double> v) : name(n), values(std::move(v)) {}

    std::string name;
    std::vector<double> values;
  };

  std::vector< NamedField > midpoint_data;

  std::string get_path(const std::string & filename, std::size_t time_idx) const;
};

class GraphData {
public:
  GraphData();

  double get_length(const Edge &) const;
  double get_A0(const Edge &) const;
  double get_G0(const Edge &) const;
  double get_radius(const Edge &) const;
  double get_C(const Edge &) const;

  double set_length(std::size_t type, double value);
  double set_A0(std::size_t type, double value);
  double set_G0(const Edge &);
  double set_radius(const Edge &);
  double set_C(const Edge &);

private:
  std::vector<std::size_t> p_id_to_vessel_type;

  std::vector<double> p_lengths;

  std::vector<double> p_A0;

  std::vector<double> p_G0;

  std::vector<double> p_radii;

  std::vector<double> p_C;
};

class GraphBoundaryData {
public:
  GraphBoundaryData();

  bool is_inflow(const Vertex &) const;

  bool is_outflow(const Vertex &) const;
};

} // namespace macrocirculation

#endif
