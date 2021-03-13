////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_STORAGE_H
#define TUMORMODELS_GRAPH_STORAGE_H

#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace macrocirculation {

class GraphStorage;

struct Point {
  Point(double x, double y, double z);

  double x, y, z;
};

Point convex_combination(const Point &left, const Point &right, double theta);

struct PhysicalData {
  PhysicalData(double G0, double A0, double rho, double length);

  // physical parameters
  double G0;
  double A0;
  double rho;
  double length;
};

struct DiscretizationData {
  std::vector<double> lengths;

  std::size_t num_micro_edges() const;
};

struct EmbeddingData {
  std::vector<Point> points;
};

class MicroPrimitive {
public:
  MicroPrimitive(std::size_t local_id, std::size_t global_id)
      : d_local_id(local_id),
        d_global_id(global_id) {}

  std::size_t get_local_id() const { return d_local_id; }
  std::size_t get_global_id() const { return d_local_id; }

private:
  std::size_t d_local_id;
  std::size_t d_global_id;
};

using MicroVertex = MicroPrimitive;
using MicroEdge = MicroPrimitive;

class Primitive {
public:
  explicit Primitive(std::size_t id);

  std::size_t get_id() const;

protected:
  std::size_t p_id;
};

class Vertex : public Primitive {
public:
  explicit Vertex(std::size_t id);

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
  bool p_inflow;

  std::function<double(double)> p_inflow_value;

  std::vector<std::size_t> p_neighbors;

  friend GraphStorage;
};

class Edge : public Primitive {
public:
  const std::vector<std::size_t> &get_vertex_neighbors() const;

  const PhysicalData &get_physical_data() const;
  const DiscretizationData &get_discretization_data() const;
  const EmbeddingData &get_embedding_data() const;

  void add_physical_data(const PhysicalData &data);
  void add_discretization_data(const DiscretizationData &data);
  void add_embedding_data(const EmbeddingData &data);

  bool has_physical_data() const;
  bool has_discretization_data() const;
  bool has_embedding_data() const;

  /*! @brief Is the edge pointing towards the given vertex? */
  bool is_pointing_to(std::size_t vertex_id) const;

  /*! @brief Returns the rank to which the dof on the given edge are assigned. */
  int rank() const;

  std::size_t num_micro_edges() const;

  const std::vector<MicroPrimitive> &micro_edges() const;

  const std::vector<MicroPrimitive> &micro_vertices() const;

  const MicroPrimitive& left_micro_vertex() const;
  const MicroPrimitive& right_micro_vertex() const;

  Edge(std::size_t id,
       const Vertex &v1,
       const Vertex &v2,
       std::size_t first_micro_edge_id,
       std::size_t first_micro_vertex_id,
       std::size_t num_micro_edges);

protected:
  /*! @brief Assigns the edge to a certain rank.
    *
    *  This method is protected, since the GraphStorage structure could cache data for each rank.
    *  Hence, the storage has to be aware of changes to invalidate its caches and therefore
    *  manages the assignment to ranks.
    */
  void assign_to_rank(int rank);

  std::vector<std::size_t> p_neighbors;

  std::unique_ptr<PhysicalData> physical_data;
  std::unique_ptr<DiscretizationData> discretization_data;
  std::unique_ptr<EmbeddingData> embedding_data;

  int d_rank;

  std::vector<MicroPrimitive> d_micro_edges;
  std::vector<MicroPrimitive> d_micro_vertices;

  friend GraphStorage;
};

class GraphStorage {
public:
  GraphStorage();

  std::shared_ptr<Edge> get_edge(std::size_t id);
  std::shared_ptr<const Edge> get_edge(std::size_t id) const;

  std::shared_ptr<Vertex> get_vertex(std::size_t id);
  std::shared_ptr<const Vertex> get_vertex(std::size_t id) const;

  std::shared_ptr<Vertex> create_vertex();
  std::shared_ptr<Edge> connect(Vertex &v1, Vertex &v2, std::size_t num_micro_edges);
  void remove(Edge &e);

  std::vector<std::size_t> get_edge_ids() const;
  std::vector<std::size_t> get_vertex_ids() const;

  double num_vertices() const;
  double num_edges() const;

  /*! @brief Returns all the edge ids assigned to the given rank. */
  std::vector<std::size_t> get_active_edge_ids(int rank) const;

  /*! @brief Returns all the edge ids for the ghost layer of main_rank w.r.t. ghost_rank.
    *
    * @param main_rank   The rank, whose ghost layer we want to know.
    * @param ghost_rank  The rank which acts as a ghost.
    * @return            A list of edge ids of the ghost.
    */
  std::vector<std::size_t> get_ghost_edge_ids(int main_rank, int ghost_rank) const;

  /*! @brief Returns all the vertex ids, which he on active neighbor edge. */
  std::vector<std::size_t> get_active_vertex_ids(int rank) const;

  void assign_edge_to_rank(Edge &edge, int rank);

private:
  std::shared_ptr<Edge> connect(Vertex &v1, Vertex &v2, std::size_t edge_id, std::size_t num_micro_edges);

  /// @brief
  ///
  /// Reorders the edges such that the first one points towards the vertex,
  /// while the second one points away from the vertex.
  /// If both edges point towards the vertex or away, an exception is thrown.
  /// If it has more than two edges, an exception is thrown.
  void reorder_edges(Vertex &v);

  /*! @brief Returns true if the given edge has a neighbor edge, which is assigned to the given rank. */
  bool edge_is_neighbor_of_rank(const Edge &e, int rank) const;

  /*! @brief Returns true if the given vertex has a neighbor edge, which is assigned to the given rank. */
  bool vertex_is_neighbor_of_rank(const Vertex &e, int rank) const;

  std::map<std::size_t, std::shared_ptr<Edge>> p_edges;
  std::map<std::size_t, std::shared_ptr<Vertex>> p_vertices;

  std::size_t p_next_edge_id;
  std::size_t p_next_vertex_id;

  std::size_t d_num_micro_edges;
  std::size_t d_num_micro_vertices;
};

} // namespace macrocirculation

#endif
