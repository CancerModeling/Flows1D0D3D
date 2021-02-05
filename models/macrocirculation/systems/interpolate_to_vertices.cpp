////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "interpolate_to_vertices.hpp"

#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {


template < std::size_t DEGREE >
void interpolate_to_vertices(const GraphStorage &graph,
                             const DofMapNetwork &map,
                             const FETypeNetwork<DEGREE> &fe,
                             const std::size_t component,
                             const std::vector<double> &dof_vector,
                             std::vector<double> &interpolated) {

  interpolated.resize(2 * graph.num_edges(), 0);

  std::vector<std::size_t> dof_indices(DEGREE+1);
  std::vector<double> dof_vector_local(DEGREE+1);
  std::vector<double> evaluated_at_qps(fe.get_phi()[0].size());

  for (auto e_id : graph.get_edge_ids()) {
    auto edge = graph.get_edge(e_id);
    map.dof_indices(*edge, dof_indices, component);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_quadrature_points(dof_vector_local, evaluated_at_qps);
    interpolated[ e_id*2 + 0] = evaluated_at_qps[0];
    interpolated[ e_id*2 + 1] = evaluated_at_qps[1];
  }
}

// template instantiations
template void interpolate_to_vertices<0>(const GraphStorage &, const DofMapNetwork &, const FETypeNetwork<0> &, const std::size_t, const std::vector<double> &, std::vector<double> &);
template void interpolate_to_vertices<1>(const GraphStorage &, const DofMapNetwork &, const FETypeNetwork<1> &, const std::size_t, const std::vector<double> &, std::vector<double> &);
template void interpolate_to_vertices<2>(const GraphStorage &, const DofMapNetwork &, const FETypeNetwork<2> &, const std::size_t, const std::vector<double> &, std::vector<double> &);
template void interpolate_to_vertices<3>(const GraphStorage &, const DofMapNetwork &, const FETypeNetwork<3> &, const std::size_t, const std::vector<double> &, std::vector<double> &);

} // namespace macrocirculation
