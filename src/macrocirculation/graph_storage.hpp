////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_STORAGE_H
#define TUMORMODELS_GRAPH_STORAGE_H

#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <string>

namespace macrocirculation {

class GraphStorage;
class Edge;
class Vertex;
class MicroEdge;
class MicroVertex;

struct Point {
  Point(double x, double y, double z);

  double x, y, z;

  static double distance(const Point &a, const Point &b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
  }
};

Point convex_combination(const Point &left, const Point &right, double theta);

struct PhysicalData {
  PhysicalData(double E, double G0, double A0, double rho, double length, double viscosity, double gamma, double radius);

  /*! @brief Initializes the physical data from a common set of "sensible" physical units.
   *
   * @param elastic_modulus   Elastic modulus of vessel in Pa.
   * @param wall_thickness    Wall-thickness of vessel in cm.
   * @param density           Blood density in kg cm^{-3}! E.g. water with 997 kg/m^3 becomes 0.997e-3 kg/cm^3.
   * @param raidus            Radius in cm.
   */
  static PhysicalData set_from_data(double elastic_modulus, double wall_thickness, double density, double gamma, double radius, double length);

  // physical parameters:

  double elastic_modulus;

  /*! @brief Prefactor for static pressure, i.e. p = G0*(sqrt{A/A0} - 1) in kg s^{-2} cm^{-1}. */
  double G0;

  /*! @brief Initial vessel area in cm. */
  double A0;

  /*! @brief Initial density in kg cm^{-3}. */
  double rho;

  /*! @brief Initial length in cm. */
  double length;

  /*! @brief Viscosity of the blood flow in [kg s^{-1} cm^{-1}]. */
  double viscosity;

  /*! @brief Shape of the bloodflow profile. */
  double gamma;

  /*! @brief Vessel radius. */
  double radius;

  double get_c0() const;
};

struct DiscretizationData {
  std::vector<double> lengths;

  std::size_t num_micro_edges() const;
};

struct EmbeddingData {
  std::vector<Point> points;
};

// TODO: Change class name to RCR-model
struct PeripheralVesselData {
  double resistance;
  double compliance;

  /*! @brief The pressure at the output of the RCR model. */
  double p_out;
};

struct VesselTreeData {
  std::vector<double> resistances;
  std::vector<double> capacitances;
  std::vector<double> radii;
  double p_out;

  /*! @brief The number of furcations at each level. 2 gives a symmetric binary tree, 1 a line. */
  size_t furcation_number;
};

struct RCLModel {
  std::vector<double> resistances;
  std::vector<double> capacitances;
  std::vector<double> inductances;

  double p_out;
};

struct LinearCharacteristicData {
  double C;
  double L;
  bool points_towards_vertex;
  double p;
  double q;
};

struct NonlinearCharacteristicData {
  double G0;
  double A0;
  double rho;
  bool points_towards_vertex;
  double p;
  double q;
};

/*! @brief Connects a vertex in one graph weakly to a vertex of a different disjoint graph. */
class InterGraphConnection {
public:
  InterGraphConnection(const std::shared_ptr<GraphStorage> &graph, const Vertex &vertex);

  const GraphStorage &get_graph() const;
  const Vertex &get_vertex() const;

private:
  std::weak_ptr<GraphStorage> d_graph;
  size_t d_vertex_id;
};

class MicroEdge {
public:
  MicroEdge(std::size_t local_id, std::size_t global_id)
      : d_local_id(local_id),
        d_global_id(global_id) {}

  std::size_t get_local_id() const { return d_local_id; }
  std::size_t get_global_id() const { return d_local_id; }

private:
  std::size_t d_local_id;
  std::size_t d_global_id;
};

class MicroVertex {
public:
  MicroVertex(std::size_t local_id, std::size_t global_id)
      : d_local_id(local_id),
        d_global_id(global_id),
        d_left_edge(nullptr),
        d_right_edge(nullptr),
        d_left_vertex(nullptr),
        d_right_vertex(nullptr) {}

  std::size_t get_local_id() const { return d_local_id; }
  std::size_t get_global_id() const { return d_local_id; }

  MicroEdge *get_left_edge() const {
    assert(d_left_edge != nullptr);
    return d_left_edge;
  }

  MicroEdge *get_right_edge() const {
    assert(d_right_edge != nullptr);
    return d_right_edge;
  }

  MicroVertex *get_left_vertex() const {
    assert(d_left_vertex != nullptr);
    return d_left_vertex;
  }
  MicroVertex *get_right_vertex() const {
    assert(d_right_vertex != nullptr);
    return d_right_vertex;
  }

protected:
  std::size_t d_local_id;
  std::size_t d_global_id;

  MicroEdge *d_left_edge;
  MicroEdge *d_right_edge;

  MicroVertex *d_left_vertex;
  MicroVertex *d_right_vertex;

  friend Edge;
};

class Primitive {
public:
  explicit Primitive(std::size_t id);

  std::size_t get_id() const;

  const std::string &get_name() const;
  void set_name(const std::string &name);

protected:
  std::size_t p_id;

  std::string p_name;
};

// TODO: add securicty assertion to Vertex for a flow type
enum class FlowType {
  Undefined,
  InflowFixedFlow,
  FreeOutflow,
  Windkessel,
  VesselTree,
  LinearCharacteristic,
  NonlinearCharacteristic,
  RCLModel,
  InflowFixedPressure
};

class Vertex : public Primitive {
public:
  explicit Vertex(std::size_t id);

  size_t local_edge_index(const Edge& ege) const;

  const std::vector<std::size_t> &get_edge_neighbors() const;

  bool is_leaf() const;

  bool is_unconnected() const;

  bool is_bifurcation() const;

  /*! @brief Returns true if the boundary conditions on this vertex have been finalized and are not allowed to change anymore. */
  bool bc_finalized() const;

  /*! @brief Finalizes the boundary conditions on this vertex. After this, they cannot be changed anymore! */
  void finalize_bcs();

  /*! @brief Marks the given vertex as part of the inflow boundary, where the given time dependent function provides the boundary values. */
  void set_to_inflow_with_fixed_flow(std::function<double(double)> inflow_value);

  /*! @brief Marks the given vertex as part of the inflow boundary, where the given time dependent function provides the boundary values. */
  void set_to_inflow_with_fixed_pressure(std::function<double(double)> pressure_value);

  /*! @brief Marks the given vertex as part of the free outflow boundary. */
  void set_to_free_outflow();

  /*! @brief Marks the given vertex as part of the free outflow boundary. */
  void set_to_windkessel_outflow(double r, double c);

  void set_to_vessel_tree_outflow(double p, const std::vector<double> &resistances, const std::vector<double> &capacitances, const std::vector<double> &radii, size_t furcation_number);

  void set_to_vessel_rcl_outflow(double p, const std::vector<double> &resistances, const std::vector<double> &capacitances, const std::vector<double> &inductances);

  void set_to_linear_characteristic_inflow(double C, double L, bool points_towards_vertex, double p, double q);

  void update_linear_characteristic_inflow(double p, double q);

  /*! @brief Updates the outgoing pressure at the vessel tips.
   *         If the boundary conditions do not have a pressure parameter a runtime_error is thrown.
   */
  void update_vessel_tip_pressures(double p);

  void set_to_nonlinear_characteristic_inflow(double G0, double A0, double rh0, bool points_towards_vertex, double p, double q);

  void update_nonlinear_characteristic_inflow(double p, double q);

  /*! @brief Returns whether the given vertex is part of the inflow boundary where the flow is given. */
  bool is_inflow_with_fixed_flow() const;

  /*! @brief Returns whether the given vertex is part of the inflow boundary where the pressure is given. */
  bool is_inflow_with_fixed_pressure() const;

  /*! @brief Returns the inflow value at the given vertex. */
  double get_inflow_value(double time) const;

  const PeripheralVesselData &get_peripheral_vessel_data() const;

  const VesselTreeData &get_vessel_tree_data() const;

  const LinearCharacteristicData &get_linear_characteristic_data() const;

  const NonlinearCharacteristicData &get_nonlinear_characteristic_data() const;

  const RCLModel &get_rcl_data() const;

  bool is_free_outflow() const;

  bool is_windkessel_outflow() const;

  bool is_vessel_tree_outflow() const;

  bool is_linear_characteristic_inflow() const;

  bool is_nonlinear_characteristic_inflow() const;

  bool is_rcl_outflow() const;

  void add_inter_graph_connection(std::shared_ptr<GraphStorage> graph, Vertex &v);

  const std::vector<InterGraphConnection> &get_inter_graph_connections() const;

  /*! @brief Establishes a symmetric inter graph connections between vertex v1 on graph g1 with vertex v2 on graph g2. */
  static void connect(const std::shared_ptr<GraphStorage> &g1, Vertex &v1, const std::shared_ptr<GraphStorage> &g2, Vertex &v2);

private:
  std::function<double(double)> p_inflow_value;

  FlowType p_flow_type;

  PeripheralVesselData p_peripheral_vessel_data;

  VesselTreeData p_vessel_tree_data;

  LinearCharacteristicData p_linear_characteristic_data;

  NonlinearCharacteristicData p_nonlinear_characteristic_data;

  RCLModel p_rcl_data;

  std::vector<InterGraphConnection> d_inter_graph_connections;

  std::vector<std::size_t> p_neighbors;

  bool d_bcs_finalized;

  friend GraphStorage;
};

class Edge : public Primitive {
public:
  struct InnerVerticesIterator {
    explicit InnerVerticesIterator(const std::vector<MicroVertex> &vertices)
        : d_vertices(vertices) {
      assert(d_vertices.size() > 2);
    }

    const MicroVertex *begin() const { return (&d_vertices.front()) + 1; }
    const MicroVertex *end() const { return (&d_vertices.back()); }

  private:
    const std::vector<MicroVertex> &d_vertices;
  };

  const std::vector<std::size_t> &get_vertex_neighbors() const;

  const PhysicalData &get_physical_data() const;
  const DiscretizationData &get_discretization_data() const;
  const EmbeddingData &get_embedding_data() const;

  void add_physical_data(const PhysicalData &data);
  void add_discretization_data(const DiscretizationData &data);
  void add_embedding_data(const EmbeddingData &data);

  std::size_t get_adajcent_micro_edge_id(const Vertex &vertex) const;

  bool has_physical_data() const;
  bool has_discretization_data() const;
  bool has_embedding_data() const;

  /*! @brief Is the edge pointing towards the given vertex? */
  bool is_pointing_to(std::size_t vertex_id) const;

  /*! @brief Returns the rank to which the dof on the given edge are assigned. */
  int rank() const;

  std::size_t num_micro_edges() const;
  std::size_t num_micro_vertices() const;

  const std::vector<MicroEdge> &micro_edges() const;

  const std::vector<MicroVertex> &micro_vertices() const;

  InnerVerticesIterator inner_micro_vertices() const;

  const MicroVertex &left_micro_vertex() const;
  const MicroVertex &right_micro_vertex() const;

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

  std::vector<MicroEdge> d_micro_edges;
  std::vector<MicroVertex> d_micro_vertices;

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

  size_t num_vertices() const;
  size_t num_edges() const;

  /*! @brief Finalizes all the boundary conditions on the vertices. */
  void finalize_bcs();

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

  std::vector<std::size_t> get_active_and_connected_vertex_ids(int rank) const;

  bool owns_primitive(const Vertex &vertex, size_t rank) const;

  void assign_edge_to_rank(Edge &edge, int rank);

  std::vector<std::shared_ptr<Vertex>> find_embedded_vertices(const Point &p) const;

  /*! @brief Searches the given named edge. */
  std::shared_ptr<Edge> find_edge_by_name(const std::string &name);

  /*! @brief Searches the given named vertex. */
  std::shared_ptr<Vertex> find_vertex_by_name(const std::string &name);

  bool has_named_vertex(const std::string &name) const;

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

  /*! @brief Returns true if the given edge has a neighbor edge on a connected graph, which is assigned to the given rank. */
  bool edge_is_connected_to_rank(const Edge &e, int rank) const;

  /*! @brief Returns true if the given vertex is to connected to another graph with an edge assigned to the given rank. */
  bool vertex_is_connected_to_rank(const Vertex &vertex, int rank) const;

  std::map<std::size_t, std::shared_ptr<Edge>> p_edges;
  std::map<std::size_t, std::shared_ptr<Vertex>> p_vertices;

  std::size_t p_next_edge_id;
  std::size_t p_next_vertex_id;

  std::size_t d_num_micro_edges;
  std::size_t d_num_micro_vertices;
};

std::vector<double> get_normals(const GraphStorage &graph, const Vertex &v);

void get_normals(const GraphStorage &graph, const Vertex &v, std::vector<double> &sigma);

} // namespace macrocirculation

#endif
