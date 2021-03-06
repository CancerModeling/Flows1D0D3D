////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "interpolate_to_vertices.hpp"

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"
#include "communication/mpi.hpp"

namespace macrocirculation {

template<std::size_t DEGREE>
void interpolate_to_vertices(const MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMapNetwork &map,
                             const std::size_t component,
                             const std::vector<double> &dof_vector,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated) {

  interpolated.resize(2 * graph.num_edges(), 0);

  FETypeNetwork<DEGREE> fe(create_trapezoidal_rule());

  std::vector<std::size_t> dof_indices(DEGREE + 1);
  std::vector<double> dof_vector_local(DEGREE + 1);
  std::vector<double> evaluated_at_qps(fe.get_phi()[0].size());

  const auto rank = mpi::rank(comm);

  for (auto e_id : graph.get_edge_ids()) {
    auto edge = graph.get_edge(e_id);

    // we only interpolate on the given graph
    if (edge->rank() != rank)
      continue;

    map.dof_indices(*edge, dof_indices, component);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_quadrature_points(dof_vector_local, evaluated_at_qps);
    interpolated[e_id * 2 + 0] = evaluated_at_qps[0];
    interpolated[e_id * 2 + 1] = evaluated_at_qps[1];

    points[e_id * 2 + 0] = edge->get_coordinate_v0();
    points[e_id * 2 + 1] = edge->get_coordinate_v1();
  }
}

// template instantiations

template void interpolate_to_vertices<0>(const MPI_Comm, const GraphStorage &, const DofMapNetwork &, const std::size_t, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void interpolate_to_vertices<1>(const MPI_Comm, const GraphStorage &, const DofMapNetwork &, const std::size_t, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void interpolate_to_vertices<2>(const MPI_Comm, const GraphStorage &, const DofMapNetwork &, const std::size_t, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void interpolate_to_vertices<3>(const MPI_Comm, const GraphStorage &, const DofMapNetwork &, const std::size_t, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);

} // namespace macrocirculation
