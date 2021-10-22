////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "quantities_of_interest.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "interpolate_to_vertices.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

void calculate_total_pressure(const MPI_Comm comm,
                              const GraphStorage &graph,
                              const DofMap &map,
                              const std::vector<double> &dof_vector,
                              std::vector<Point> &points,
                              std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    const auto local_dof_map = map.get_local_dof_map(*edge);

    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions());
    std::vector<double> dof_vector_local(local_dof_map.num_basis_functions());

    FETypeNetwork fe(create_trapezoidal_rule(), local_dof_map.num_basis_functions() - 1);

    if (embedding.points.size() == 2 && local_dof_map.num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], local_dof_map.num_micro_edges(), points);
    else if (embedding.points.size() == local_dof_map.num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      fe.reinit(h);

      local_dof_map.dof_indices(micro_edge_id, 0, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto Q = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      local_dof_map.dof_indices(micro_edge_id, 1, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto A = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      interpolated.push_back(nonlinear::get_p_from_QA(Q.left, A.left, param));
      interpolated.push_back(nonlinear::get_p_from_QA(Q.right, A.right, param));
    }
  }
}

void calculate_static_pressure(const MPI_Comm comm,
                               const GraphStorage &graph,
                               const DofMap &map,
                               const std::vector<double> &dof_vector,
                               std::vector<Point> &points,
                               std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    const auto local_dof_map = map.get_local_dof_map(*edge);

    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions());
    std::vector<double> dof_vector_local(local_dof_map.num_basis_functions());

    FETypeNetwork fe(create_trapezoidal_rule(), local_dof_map.num_basis_functions() - 1);

    if (embedding.points.size() == 2 && local_dof_map.num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], local_dof_map.num_micro_edges(), points);
    else if (embedding.points.size() == local_dof_map.num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      fe.reinit(h);

      local_dof_map.dof_indices(micro_edge_id, 1, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto A = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      interpolated.push_back(nonlinear::get_p_from_A(A.left, param.G0, param.A0));
      interpolated.push_back(nonlinear::get_p_from_A(A.right, param.G0, param.A0));
    }
  }
}

} // namespace macrocirculation
