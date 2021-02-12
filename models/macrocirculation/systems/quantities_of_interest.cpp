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
                              const FETypeNetwork<DEGREE> &fe,
                              const std::vector<double> &dof_vector,
                              std::vector<double> &interpolated)
{
  interpolated.resize(2 * graph.num_edges(), 0);

  std::vector<std::size_t> dof_indices(DEGREE+1);
  std::vector<double> dof_vector_local(DEGREE+1);
  std::vector< double > values(2, 0);

  for (auto e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    map.dof_indices(*edge, dof_indices, 0);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_quadrature_points(dof_vector_local, values);
    const double Q_l = values[0];
    const double Q_r = values[1];
    map.dof_indices(*edge, dof_indices, 1);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_quadrature_points(dof_vector_local, values);
    const double A_l = values[0];
    const double A_r = values[1];
    const auto& data = vessel_data.get_parameters(*edge);
    interpolated[ e_id*2 + 0] = calculate_p_from_QA(Q_l, A_l, data.G0, data.rho, data.A0) ;
    interpolated[ e_id*2 + 1] = calculate_p_from_QA(Q_r, A_r, data.G0, data.rho, data.A0) ;
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
  std::vector< double > values(2, 0);

  for (auto e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);
    map.dof_indices(*edge, dof_indices, 0);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_boundary_points(dof_vector_local, values);
    const auto& data = vessel_data.get_parameters(*edge);
    interpolated[ e_id*2 + 0] = calculate_static_p(values[0], data.G0, data.A0) ;
    interpolated[ e_id*2 + 1] = calculate_static_p(values[1], data.G0, data.A0) ;
  }
}

template void calculate_total_pressure<0>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeNetwork<0> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<1>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeNetwork<1> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<2>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeNetwork<2> &, const std::vector<double> &, std::vector<double> &);
template void calculate_total_pressure<3>(const GraphStorage &, const VesselDataStorage& , const DofMapNetwork &, const FETypeNetwork<3> &, const std::vector<double> &, std::vector<double> &);

}

