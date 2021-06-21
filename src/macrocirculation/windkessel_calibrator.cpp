////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>

#include "communication/mpi.hpp"
#include "graph_storage.hpp"
#include "vessel_formulas.hpp"
#include "windkessel_calibrator.hpp"

namespace macrocirculation {

WindkesselCalibrator::WindkesselCalibrator(std::shared_ptr<GraphStorage> graph, bool verbose)
    : d_graph(std::move(graph)),
      d_verbose(verbose),
      d_total_C_edge(0),
      d_total_R(1.34),   // 1.34e8 Pa s / m^3   ([Pa s / m^-3] = 10^{-8} [kg/s cm])
      d_total_C(9.45e-1) // 9.e5e-9 m^3/Pa ([m^3 / Pa] = 10^8 [cm s^2 / kg])
{
  reset();
}

void WindkesselCalibrator::reset() {
  reset_total_flows();
  calculate_total_edge_capacitance();
}

void WindkesselCalibrator::reset_total_flows() {
  d_total_flows.clear();

  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    if (d_graph->get_vertex(v_id)->is_leaf())
      d_total_flows[v_id] = 0.;
  }
}

void WindkesselCalibrator::calculate_total_edge_capacitance() {
  d_total_C_edge = 0;
  for (auto e_id : d_graph->get_active_edge_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto e = d_graph->get_edge(e_id);
    auto &data = e->get_physical_data();

    double c0 = calculate_c0(data.G0, data.rho, data.A0); // [cm/s]

    // in meter
    double C_e = (1e-2 * data.length) * (1e-4 * data.A0) / (1060 * std::pow(1e-2 * c0, 2));

    d_total_C_edge += C_e;
  }

  MPI_Allreduce(&d_total_C_edge, &d_total_C_edge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

std::map<size_t, RCRData> WindkesselCalibrator::estimate_parameters_local() {

  // print the total flows
  double sum_of_flows = 0;
  double sum_of_out_flows = 0;
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_leaf()) {
      if (d_verbose)
        std::cout << "flow at " << v_id << " = " << d_total_flows[v_id] << std::endl;
      sum_of_flows += d_total_flows[v_id];
      if (v->is_free_outflow()) {
        sum_of_out_flows += d_total_flows[v_id];
      }
    }
  }
  MPI_Allreduce(&sum_of_flows, &sum_of_flows, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sum_of_out_flows, &sum_of_out_flows, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (d_verbose)
    std::cout << "flow-sum = " << sum_of_flows << std::endl;

  std::map<size_t, RCRData> resistances;

  // divide resistances
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_free_outflow()) {
      auto z = d_total_flows[v_id] / sum_of_out_flows;
      double local_R = d_total_R / z;
      double C = z * (d_total_C - d_total_C_edge);

      if (d_verbose)
        std::cout << "vertex=" << v->get_name() << " R=" << local_R * 1e-1 << "10^9 Pa s m^{-3}, C=" << C * 1e2 << " 10^{-10} m^3 Pa^{-1}" << std::endl;

      resistances[v_id] = {local_R, C};
    }
  }

  return resistances;
};

std::map<size_t, RCRData> WindkesselCalibrator::estimate_parameters() {
  auto local_parameters = estimate_parameters_local();
  std::map<size_t, RCRData> global_parameters;

  std::array<double, 2> data{0., 0.};

  for (auto v_id : d_graph->get_vertex_ids()) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_free_outflow()) {
      auto &e = *d_graph->get_edge(v->get_edge_neighbors()[0]);
      if (e.rank() == mpi::rank(MPI_COMM_WORLD)) {
        data[0] = local_parameters[v_id].resistance;
        data[1] = local_parameters[v_id].capacitance;
      }
      MPI_Bcast(&data[0], 2, MPI_DOUBLE, e.rank(), MPI_COMM_WORLD);
      global_parameters[v_id] = {data[0], data[1]};
    }
  }

  return global_parameters;
}

void WindkesselCalibrator::set_total_C(double C) { d_total_C = C; }

void WindkesselCalibrator::set_total_R(double R) { d_total_R = R; }

void parameters_to_json(const std::string &filepath, const std::map<size_t, RCRData> &parameters, const std::shared_ptr<GraphStorage> &storage) {
  // only root writes the file
  if (mpi::rank(MPI_COMM_WORLD) != 0)
    return;

  using json = nlohmann::json;

  json j;

  j["vertices"] = json::array();

  for (auto it : parameters) {
    const auto v_id = it.first;
    const auto data = it.second;
    const auto &v = *storage->get_vertex(v_id);

    if (v.get_name().empty())
      throw std::runtime_error("cannot serialize unnamed vertex");

    j["vertices"].push_back({{"name", v.get_name()},
                             {"peripheral_resistance", data.resistance},
                             {"peripheral_compliance", data.capacitance}});
  }

  // write back
  std::ofstream o(filepath);
  o << std::setw(4);
  o << j;
}


} // namespace macrocirculation
