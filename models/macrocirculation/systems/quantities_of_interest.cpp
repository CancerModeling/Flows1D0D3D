////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "quantities_of_interest.hpp"

#include "graph_storage.hpp"
#include "vessel_data_storage.hpp"
#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

template < std::size_t DEGREE >
void calculate_total_pressure(const GraphStorage &graph,
                              const VesselDataStorage& vessel_data,
                              const DofMapNetwork &map,
                              const FETypeInnerBdryNetwork<DEGREE> &fe,
                              const std::vector<double> &dof_vector,
                              std::vector<double> &interpolated)
{
  interpolated.resize(2 * graph.num_edges(), 0);

  std::vector<std::size_t> dof_indices(DEGREE+1);
  std::vector<double> dof_vector_local(DEGREE+1);

  for (auto e_id : graph.get_edge_ids()) {
    auto edge = graph.get_edge(e_id);
    map.dof_indices(*edge, dof_indices, 0);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    auto Q_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);
    map.dof_indices(*edge, dof_indices, 1);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    auto A_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);
    auto data = vessel_data.get_parameters(*edge);
    interpolated[ e_id*2 + 0] = calculate_p_from_QA(Q_values.left, A_values.left, data.G0, data.rho, data.A0) ;
    interpolated[ e_id*2 + 1] = calculate_p_from_QA(Q_values.right, A_values.right, data.G0, data.rho, data.A0) ;
  }
}

template < std::size_t DEGREE >
void calculate_static_pressure(const GraphStorage &graph,
                              const VesselDataStorage& vessel_data,
                              const DofMapNetwork &map,
                              const FETypeInnerBdryNetwork<DEGREE> &fe,
                              const std::vector<double> &dof_vector,
                              std::vector<double> &interpolated)
{
  interpolated.resize(2 * graph.num_edges(), 0);

  std::vector<std::size_t> dof_indices(DEGREE+1);
  std::vector<double> dof_vector_local(DEGREE+1);

  for (auto e_id : graph.get_edge_ids()) {
    auto edge = graph.get_edge(e_id);
    map.dof_indices(*edge, dof_indices, 0);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    auto A_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);
    auto data = vessel_data.get_parameters(*edge);
    interpolated[ e_id*2 + 0] = calculate_static_p(A_values.left, data.G0, data.A0) ;
    interpolated[ e_id*2 + 1] = calculate_static_p(A_values.right, data.G0, data.A0) ;
  }
}

template void calculate_total_pressure<0>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeInnerBdryNetwork<0> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<1>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeInnerBdryNetwork<1> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<2>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeInnerBdryNetwork<2> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<3>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeInnerBdryNetwork<3> &, const std::vector<double> &, std::vector<double> &);

}

