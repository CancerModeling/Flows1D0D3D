////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_STORAGE_H
#define TUMORMODELS_GRAPH_STORAGE_H

#include <functional>
#include <libmesh/point.h>
#include <map>
#include <memory>
#include <utility>
#include <vector>

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
  bool is_bifurcation() const;

  /*! @brief Marks the given vertex as part of the inflow boundary, where the given time dependent function provides the boundary values. */
  void set_to_inflow(std::function<double(double)> inflow_value);

  /*! @brief Marks the given vertex as part of the outflow boundary. */
  void set_to_outflow();

  /*! @brief Returns whether the given vertex is part of the inflow boundary. */
  bool is_inflow() const;

  /*! @brief Returns the inflow value at the given vertex. */
  double get_inflow_value(double time) const;

private:
  Point p_coordinate;

  bool p_inflow;

  std::function<double(double)> p_inflow_value;

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

  /*! @brief Is the edge pointing towards the given vertex? */
  bool is_pointing_to(std::size_t vertex_id) const;

  std::size_t get_type_id() const;

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

  double num_vertices() const;
  double num_edges() const;

  /*! @brief Adds a line segment starting at "from" and ending at "to" to the graph. */
  void line_to(Vertex &from, Vertex &to, std::size_t vessel_id, std::size_t num_edges);

private:
  std::shared_ptr<Edge> connect(Vertex &v1, Vertex &v2, std::size_t edge_type_id, std::size_t edge_id);

  /// @brief
  ///
  /// Reorders the edges such that the first one points towards the vertex,
  /// while the second one points away from the vertex.
  /// If both edges point towards the vertex or away, an exception is thrown.
  /// If it has more than two edges, an exception is thrown.
  void reorder_edges(Vertex &v);

  std::map<std::size_t, std::shared_ptr<Edge>> p_edges;
  std::map<std::size_t, std::shared_ptr<Vertex>> p_vertices;

  std::size_t p_next_edge_id;
  std::size_t p_next_vertex_id;
};

} // namespace macrocirculation

#endif
