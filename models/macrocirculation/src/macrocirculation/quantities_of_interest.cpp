////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "quantities_of_interest.hpp"

#include "communication/mpi.hpp"
#include "dof_map_network.hpp"
#include "fe_type_network.hpp"
#include "graph_storage.hpp"
#include "vessel_data_storage.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

template<std::size_t DEGREE>
void calculate_total_pressure(const MPI_Comm comm,
                              const GraphStorage &graph,
                              const VesselDataStorage &vessel_data,
                              const DofMapNetwork &map,
                              const std::vector<double> &dof_vector,
                              std::vector<Point> &points,
                              std::vector<double> &interpolated) {
  if (points.size() != interpolated.size() || interpolated.size() != 2 * graph.num_edges())
    throw std::runtime_error("vectors have wrong size.");

  std::vector<std::size_t> dof_indices(DEGREE + 1);
  std::vector<double> dof_vector_local(DEGREE + 1);
  std::vector<double> values(2, 0);

  FETypeNetwork<DEGREE> fe(create_trapezoidal_rule());

  const auto rank = mpi::rank(comm);

  for (auto e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);

    // we only calculate on the given graph
    if (edge->rank() != rank)
      continue;

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
    const auto &data = vessel_data.get_parameters(*edge);
    interpolated[e_id * 2 + 0] = calculate_p_from_QA(Q_l, A_l, data.G0, data.rho, data.A0);
    interpolated[e_id * 2 + 1] = calculate_p_from_QA(Q_r, A_r, data.G0, data.rho, data.A0);
    points[e_id * 2 + 0] = edge->get_coordinate_v0();
    points[e_id * 2 + 1] = edge->get_coordinate_v1();
  }
}

template<std::size_t DEGREE>
void calculate_static_pressure(const MPI_Comm comm,
                               const GraphStorage &graph,
                               const VesselDataStorage &vessel_data,
                               const DofMapNetwork &map,
                               const std::vector<double> &dof_vector,
                               std::vector<Point> &points,
                               std::vector<double> &interpolated) {
  if (points.size() != interpolated.size() || interpolated.size() != 2 * graph.num_edges())
    throw std::runtime_error("vectors have wrong size.");

  std::vector<std::size_t> dof_indices(DEGREE + 1);
  std::vector<double> dof_vector_local(DEGREE + 1);
  std::vector<double> values(2, 0);

  FETypeNetwork<DEGREE> fe(create_trapezoidal_rule());

  const auto rank = mpi::rank(comm);

  for (auto e_id : graph.get_edge_ids()) {
    const auto edge = graph.get_edge(e_id);

    // we only calculate on the given graph
    if (edge->rank() != rank)
      continue;

    map.dof_indices(*edge, dof_indices, 1);
    extract_dof(dof_indices, dof_vector, dof_vector_local);
    fe.evaluate_dof_at_quadrature_points(dof_vector_local, values);
    const auto &data = vessel_data.get_parameters(*edge);
    interpolated[e_id * 2 + 0] = calculate_static_p(values[0], data.G0, data.A0);
    interpolated[e_id * 2 + 1] = calculate_static_p(values[1], data.G0, data.A0);
    points[e_id * 2 + 0] = edge->get_coordinate_v0();
    points[e_id * 2 + 1] = edge->get_coordinate_v1();
  }
}

template void calculate_total_pressure<0>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_total_pressure<1>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_total_pressure<2>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_total_pressure<3>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);

template void calculate_static_pressure<0>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_static_pressure<1>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_static_pressure<2>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);
template void calculate_static_pressure<3>(const MPI_Comm, const GraphStorage &, const VesselDataStorage &, const DofMapNetwork &, const std::vector<double> &, std::vector<Point> &, std::vector<double> &);

} // namespace macrocirculation
