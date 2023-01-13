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

EdgeTreeParameters calculate_edge_tree_parameters(const Edge &edge) {
  auto &param = edge.get_physical_data();
  const double E = param.elastic_modulus * 4;
  const double r_0 = param.radius;
  const double r_cap = 7.5e-4;
  const double h_0 = 0.003;
  const double gamma = 2.5;
  // floor ?
  const int N = static_cast<int>(std::ceil(gamma * std::log(r_0 / r_cap) / std::log(2)));
  const auto alpha = 1. / std::pow(2, 1 / gamma);
  std::vector<double> list_C;
  std::vector<double> list_R;
  std::vector<double> list_radii;
  std::vector<double> list_lengths;
  double r = r_0 * alpha;
  for (size_t k = 0; k < static_cast< size_t > ( N ); k += 1) {
    const double l = 280 * r;
    const double C = 3 * std::pow(r, 3) * M_PI * l / (2 * E * h_0);
    const double viscosity = viscosity_bloodplasma(r);
    double R = 2 * (param.gamma + 2) * viscosity * l / (std::pow(r, 2));
    list_C.push_back(C);
    list_R.push_back(R);
    list_radii.push_back(r);
    list_lengths.push_back(l);
    r *= alpha;
  }
  return {list_lengths, list_radii, list_R, list_C};
}

void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph) {
  set_0d_tree_boundary_conditions(graph, [](auto) { return true; });
}

void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph, const std::string &prefix) {
  set_0d_tree_boundary_conditions(graph, [&](const Vertex &v) { return v.get_name().rfind(prefix) == 0; });
}

void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph, const std::function<bool(const Vertex &)> &conditional) {
  for (auto &v_id : graph->get_vertex_ids()) {
    auto &vertex = *graph->get_vertex(v_id);
    if (!vertex.is_leaf())
      continue;

    if (vertex.is_inflow_with_fixed_flow() || vertex.is_inflow_with_fixed_pressure()) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " found inflow " << vertex.get_name() << std::endl;
      continue;
    }

    if (vertex.is_nonlinear_characteristic_inflow() || vertex.is_linear_characteristic_inflow()) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " found characteristic inflow " << vertex.get_name() << std::endl;
      continue;
    }

    // we do not touch the circle of willis
    if (!conditional(vertex)) {
      std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " ignoring node " << vertex.get_name() << " and keeps windkessel bc." << std::endl;
      continue;
    }


    std::cout << "rank = " << mpi::rank(MPI_COMM_WORLD) << " sets " << vertex.get_name() << " to tree bc" << std::endl;
    auto &edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &param = edge.get_physical_data();
    // const double E = param.elastic_modulus;
    const double E = param.elastic_modulus * 4;
    const double r_0 = param.radius;
    //const double r_cap = 7.5e-4;
    const double r_cap = 0.001;
    //const double h_0 = 1e-4;
    const double h_0 = 0.003;
    const double p_cap = 30 * (133.333) * 1e-2;
    const double gamma = 2.5;
    // floor ?
    const int N = static_cast<int>(std::ceil(gamma * std::log(r_0 / r_cap) / std::log(2)));
    const auto alpha = 1. / std::pow(2, 1 / gamma);
    // const auto beta = 3. / 4.;
    std::vector<double> list_radius;
    std::vector<double> list_C;
    std::vector<double> list_R;
    double r = r_0 * alpha;
    //double l = 0.2; // 2mm is the average vessel length
    //double l = 3 * beta; // 2mm is the average vessel length
    std::cout << "vessel start" << std::endl;
    for (size_t k = 0; k < static_cast< size_t > ( N ); k += 1) {
      const double l = 280 * r;
      std::cout << " l = " << l << " r = " << r << std::endl;
      const double C = 3 * std::pow(r, 3) * M_PI * l / (2 * E * h_0);
      const double viscosity = viscosity_bloodplasma(r);
      double R = 2 * (param.gamma + 2) * viscosity * l / (std::pow(r, 2));
      list_C.push_back(C);
      list_R.push_back(R);
      list_radius.push_back(r);
      r *= alpha;
    }
    std::cout << "vessel stop" << std::endl;

    vertex.set_to_vessel_tree_outflow(p_cap, list_R, list_C, list_radius, 2);
  }
}

void convert_rcr_to_partitioned_tree_bcs(const std::shared_ptr<GraphStorage> &graph) {
  for (auto &v_id : graph->get_vertex_ids()) {
    auto &vertex = *graph->get_vertex(v_id);

    if (vertex.is_windkessel_outflow()) {
      auto &edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &data = vertex.get_peripheral_vessel_data();

      std::cout << "rank " << mpi::rank(MPI_COMM_WORLD) << " sets vertex " << vertex.get_name()
                << " (id = " << vertex.get_id() << ", neighbor-edge-id = " << edge.get_id() << ")" << std::endl;

      const double R1 = calculate_R1(edge.get_physical_data());
      const double R2 = data.resistance - R1;
      const double C = data.compliance;

      const double p_cap = 5 * (133.333) * 1e-2;

      std::vector<double> list_C;
      std::vector<double> list_R;

      // arterioles
      list_C.push_back(C * 0.08);
      list_C.push_back(C * 0.09);
      list_C.push_back(C * 0.10);
      list_C.push_back(C * 0.10);
      list_C.push_back(C * 0.11);
      list_C.push_back(C * 0.12);
      list_C.push_back(C * 0.48);
      list_C.push_back(C * 0.68);
      // capillaries
      list_C.push_back(C * 0.14);
      // venules
      list_C.push_back(C * 0.20);
      // small veins
      list_C.push_back(C * 0.40);

      // arterioles 60
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.20);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.05);
      // capillaries
      list_R.push_back(R2 * 0.05);
      // venules
      list_R.push_back(R2 * 0.025);
      // small veins
      list_R.push_back(R2 * 0.025);

      std::vector<double> radii;
      for (int i = 0; i < static_cast< int > ( list_R.size() ); i += 1)
        radii.push_back(NAN);

      vertex.set_to_vessel_tree_outflow(p_cap, list_R, list_C, radii, 1);
    }
  }
}

void convert_rcr_to_rcl_chain_bcs(const std::shared_ptr<GraphStorage> &graph) {
  for (auto &v_id : graph->get_vertex_ids()) {
    auto &vertex = *graph->get_vertex(v_id);

    if (vertex.is_windkessel_outflow()) {
      auto &edge = *graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &data = vertex.get_peripheral_vessel_data();

      std::cout << "rank " << mpi::rank(MPI_COMM_WORLD) << " sets vertex " << vertex.get_name()
                << " (id = " << vertex.get_id() << ", neighbor-edge-id = " << edge.get_id() << ")" << std::endl;

      const double R1 = calculate_R1(edge.get_physical_data());
      const double R2 = data.resistance - R1;
      const double C = data.compliance;

      const double p_cap = 5 * (133.333) * 1e-2;

      std::vector<double> list_C;
      std::vector<double> list_R;
      std::vector<double> list_L;

      // arterioles
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      list_C.push_back(C * 0.03125);
      // capillaries
      list_C.push_back(C * 0.15);
      // venules
      list_C.push_back(C * 0.20);
      // small veins
      list_C.push_back(C * 0.40);

      // arterioles 60
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.20);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.15);
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.10);
      list_R.push_back(R2 * 0.05);
      // capillaries
      list_R.push_back(R2 * 0.05);
      // venules
      list_R.push_back(R2 * 0.025);
      // small veins
      list_R.push_back(R2 * 0.025);

      // arterioles 60
      const double L_art = 1e-8 * 1.333;
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      list_L.push_back(L_art * 0.125);
      // capillaries
      const double L_cap = 1e-8 * 1.333;
      list_L.push_back(L_cap);
      // venules
      const double L_ven = 5e-8 * 1.333;
      list_L.push_back(L_ven);
      // small veins
      list_L.push_back(L_ven);
      std::cout << " L " << L_ven << std::endl;

      vertex.set_to_vessel_rcl_outflow(p_cap, list_R, list_C, list_L);
    }
  }
}

} // namespace macrocirculation
