////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>

#include "0d_boundary_conditions.hpp"
#include "communication/mpi.hpp"
#include "graph_storage.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph, const std::string &name_prefix) {
  for (auto &v_id : graph->get_vertex_ids()) {
    auto &vertex = *graph->get_vertex(v_id);
    if (!vertex.is_leaf())
      continue;

    if (vertex.is_inflow()) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " found inflow " << vertex.get_name() << std::endl;
      continue;
    }

    if (vertex.is_nonlinear_characteristic_inflow() || vertex.is_linear_characteristic_inflow()) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " found characteristic inflow " << vertex.get_name() << std::endl;
      continue;
    }

    // we do not touch the circle of willis
    if (vertex.get_name().rfind(name_prefix) != 0) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " ignoring node " << vertex.get_name() << " and keeps windkessel bc." << std::endl;
      continue;
    }

    std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " sets " << vertex.get_name() << " to tree bc" << std::endl;
    auto &edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &param = edge.get_physical_data();
    const double E = param.elastic_modulus;
    const double r_0 = param.radius;
    const double r_cap = 5e-4;
    const double h_0 = 1e-4;
    const double p_cap = 30 * (133.333) * 1e-2;
    const int N = static_cast<int>(std::ceil(3 * std::log(r_0 / r_cap) / std::log(2.)));
    const auto alpha = 1. / std::pow(2, 1 / 3.);
    std::vector<double> list_C;
    std::vector<double> list_R;
    double r = r_0 * alpha;
    double l = 0.2; // 2mm is the average vessel length
    for (size_t k = 0; k < N; k += 1) {
      const double C = 3 * std::pow(r, 3) * M_PI * l / (2 * E * h_0);
      const double R = 2 * (param.gamma + 2) * param.viscosity * l / (std::pow(r, 2));
      list_C.push_back(C);
      list_R.push_back(R);
      r *= alpha;
    }

    vertex.set_to_vessel_tree_outflow(p_cap, list_R, list_C, 2);
  }
}

void convert_rcr_to_partitioned_tree_bcs(const std::shared_ptr<GraphStorage> &graph)
{
  for (auto &v_id : graph->get_vertex_ids()) {
    auto &vertex = *graph->get_vertex(v_id);

    if (vertex.is_windkessel_outflow())
    {
      auto& edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto& data = vertex.get_peripheral_vessel_data();

      std::cout << "rank " << mpi::rank(MPI_COMM_WORLD) << " sets vertex " << vertex.get_name()
                << " (id = " << vertex.get_id() << ", neighbor-edge-id = " << edge.get_id() << ")" << std::endl;

      const double R1 = calculate_R1(edge.get_physical_data());
      const double R2 = data.resistance - R1;
      const double C = data.compliance;

      const double p_cap = 5 * (133.333) * 1e-2;

      std::vector<double> list_C;
      std::vector<double> list_R;

      // arterioles
      list_C.push_back(C*0.10);
      list_C.push_back(C*0.07);
      list_C.push_back(C*0.05);
      list_C.push_back(C*0.03);
      // capillaries
      list_C.push_back(C*0.15);
      // venules
      list_C.push_back(C*0.20);
      // small veins
      list_C.push_back(C*0.40);

      // arterioles 60
      list_R.push_back(R2*0.10);
      list_R.push_back(R2*0.10);
      list_R.push_back(R2*0.20);
      list_R.push_back(R2*0.25);
      // capillaries
      list_R.push_back(R2*0.15);
      // venules
      list_R.push_back(R2*0.10);
      // small veins
      list_R.push_back(R2*0.05);

      vertex.set_to_vessel_tree_outflow(p_cap, list_R, list_C, 1);
    }
  }
}

} // namespace macrocirculation