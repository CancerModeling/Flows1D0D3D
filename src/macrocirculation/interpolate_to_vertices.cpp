////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "interpolate_to_vertices.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "implicit_transport_solver.hpp"
#include "petsc/petsc_vec.hpp"

namespace macrocirculation {

void linear_interpolate_points(const Point &left,
                               const Point &right,
                               const std::size_t num_micro_edges,
                               std::vector<Point> &points) {
  for (std::size_t micro_edge_id = 0; micro_edge_id < num_micro_edges; micro_edge_id += 1) {
    points.push_back(convex_combination(left, right, static_cast<double>(micro_edge_id) / num_micro_edges));
    points.push_back(convex_combination(left, right, static_cast<double>(micro_edge_id + 1) / num_micro_edges));
  }
}

void add_discontinuous_points(const std::vector<Point> &embedded_points, std::vector<Point> &points) {
  if (embedded_points.size() < 2)
    throw std::runtime_error("not enough points in embedding");

  for (std::size_t micro_edge_id = 0; micro_edge_id < embedded_points.size() - 1; micro_edge_id += 1) {
    points.push_back(embedded_points[micro_edge_id]);
    points.push_back(embedded_points[micro_edge_id + 1]);
  }
}

void interpolate_transformation(const MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMap &map,
                             const std::size_t component,
                             const PetscVec &dof_vector,
                             std::function< double(double, const Edge&) > trafo,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    auto local_dof_map = map.get_local_dof_map(*edge);

    if (embedding.points.size() == 2 && local_dof_map.num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], local_dof_map.num_micro_edges(), points);
    else if (embedding.points.size() == local_dof_map.num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    FETypeNetwork fe(create_trapezoidal_rule(), local_dof_map.num_basis_functions() - 1);

    dof_indices.resize(local_dof_map.num_basis_functions());
    dof_vector_local.resize(local_dof_map.num_basis_functions());
    evaluated_at_qps.resize(fe.num_quad_points());

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    fe.reinit(h);

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      interpolated.push_back(trafo(boundary_values.left, *edge));
      interpolated.push_back(trafo(boundary_values.right, *edge));
    }
  }
}


void interpolate_to_vertices(const MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMap &map,
                             const std::size_t component,
                             const std::vector<double> &dof_vector,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    auto local_dof_map = map.get_local_dof_map(*edge);

    if (embedding.points.size() == 2 && local_dof_map.num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], local_dof_map.num_micro_edges(), points);
    else if (embedding.points.size() == local_dof_map.num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    FETypeNetwork fe(create_trapezoidal_rule(), local_dof_map.num_basis_functions() - 1);

    dof_indices.resize(local_dof_map.num_basis_functions());
    dof_vector_local.resize(local_dof_map.num_basis_functions());
    evaluated_at_qps.resize(fe.num_quad_points());

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    fe.reinit(h);

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      interpolated.push_back(boundary_values.left);
      interpolated.push_back(boundary_values.right);
    }
  }
}

void interpolate_to_vertices(MPI_Comm comm,
                             const GraphStorage &graph,
                             const UpwindProvider &upwind,
                             double t,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    if (embedding.points.size() == 2 && edge->num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], edge->num_micro_edges(), points);
    else if (embedding.points.size() == edge->num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    std::vector< double > upwinded_values_v( edge->num_micro_vertices() );

    upwind.get_upwinded_values(t, *edge, upwinded_values_v );

    interpolated.push_back(upwinded_values_v.front());

    for (size_t k = 1; k <upwinded_values_v.size()-1; k+=1)
    {
      interpolated.push_back(upwinded_values_v[k]);
      interpolated.push_back(upwinded_values_v[k]);
    }

    interpolated.push_back(upwinded_values_v.back());

  }
}

void fill_with_radius(const MPI_Comm comm, const GraphStorage &graph, std::vector<Point> &points, std::vector<double> &interpolated) {
  auto f = [](const Edge &e) {
    if (!e.has_physical_data())
      throw std::runtime_error("cannot get radius for edges without physical parameters");
    return e.get_physical_data().radius;
  };
  fill_with_edge_parameter(comm, graph, f, points, interpolated);
}

void fill_with_vessel_id(const MPI_Comm comm,
                         const GraphStorage &graph,
                         std::vector<Point> &points,
                         std::vector<double> &interpolated) {
  auto f = [](const Edge &e) { return e.get_id(); };
  fill_with_edge_parameter(comm, graph, f, points, interpolated);
}

void fill_with_vessel_A0(const MPI_Comm comm,
                         const GraphStorage &graph,
                         std::vector<Point> &points,
                         std::vector<double> &interpolated) {
  auto f = [](const Edge &e) { return e.get_physical_data().A0; };
  fill_with_edge_parameter(comm, graph, f, points, interpolated);
}

void fill_with_vertex_id(const MPI_Comm comm,
                         const GraphStorage &graph,
                         std::vector<Point> &points,
                         std::vector<double> &interpolated) {
  auto f = [](const Edge &e, size_t micro_edge_id)-> double {
    return 2*micro_edge_id < e.num_micro_edges() ? e.get_vertex_neighbors()[0] : e.get_vertex_neighbors()[1];
  };
  fill_with_edge_parameter(comm, graph, f, points, interpolated);
}

void fill_with_edge_parameter(const MPI_Comm comm, const GraphStorage &graph, std::function<double(const Edge &)> extractor, std::vector<Point> &points, std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    if (embedding.points.size() == 2 && edge->num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], edge->num_micro_edges(), points);
    else if (embedding.points.size() == edge->num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    const double quantity = extractor(*edge);

    for (std::size_t micro_edge_id = 0; micro_edge_id < edge->num_micro_edges(); micro_edge_id += 1) {
      interpolated.push_back(quantity);
      interpolated.push_back(quantity);
    }
  }
}

void fill_with_edge_parameter(const MPI_Comm comm, const GraphStorage &graph, const std::function<double(const Edge &, size_t micro_edge_id)>& extractor, std::vector<Point> &points, std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    if (embedding.points.size() == 2 && edge->num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], edge->num_micro_edges(), points);
    else if (embedding.points.size() == edge->num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    for (std::size_t micro_edge_id = 0; micro_edge_id < edge->num_micro_edges(); micro_edge_id += 1) {
      const double quantity = extractor(*edge, micro_edge_id);
      interpolated.push_back(quantity);
      interpolated.push_back(quantity);
    }
  }
}

} // namespace macrocirculation
