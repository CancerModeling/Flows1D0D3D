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

#include "windkessel_calibrator.hpp"

#include "communication/mpi.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "graph_storage.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

FlowIntegrator::FlowIntegrator(std::shared_ptr<GraphStorage> graph)
    : d_graph(std::move(graph)) {
  reset();
}

void FlowIntegrator::reset() {
  d_total_flows.clear();

  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    if (d_graph->get_vertex(v_id)->is_free_outflow())
      d_total_flows[v_id] = 0.;
  }
}

void FlowIntegrator::update_flow(const ExplicitNonlinearFlowSolver &solver, double tau) {
  update_flow_abstract(solver, tau);
}

void FlowIntegrator::update_flow(const ImplicitLinearFlowSolver &solver, double tau) {
  update_flow_abstract(solver, tau);
}

template<typename Solver>
void FlowIntegrator::update_flow_abstract(const Solver &solver, double tau) {
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_leaf()) {
      auto& edge = *d_graph->get_edge(v->get_edge_neighbors()[0]);
      double sigma = edge.is_pointing_to(v_id) ? +1. : -1;
      double p, q;
      solver.get_1d_pq_values_at_vertex(*v, p, q);
      d_total_flows[v_id] += sigma * q * tau;
    }
  }
}

FlowData FlowIntegrator::get_free_outflow_data() const {
  FlowData data;
  data.total_flow = 0;

  for (auto v_id : d_graph->get_vertex_ids()) {
    auto v = d_graph->get_vertex(v_id);
    if (v->is_free_outflow()) {
      auto &e = *d_graph->get_edge(v->get_edge_neighbors()[0]);
      double q = 0;
      if (e.rank() == mpi::rank(MPI_COMM_WORLD)) {
        q = d_total_flows.at(v_id);
      }
      std::cout << v_id << " " << q << std::endl;
      MPI_Bcast(&q, 1, MPI_DOUBLE, e.rank(), MPI_COMM_WORLD);
      data.flows[v_id] = q;
      data.total_flow += q;
    }
  }
  return data;
};

double get_total_edge_capacitance(const std::shared_ptr<GraphStorage> &graph) {
  double total_C_edge = 0;
  for (auto e_id : graph->get_active_edge_ids(mpi::rank(MPI_COMM_WORLD))) {
    auto e = graph->get_edge(e_id);
    auto &data = e->get_physical_data();

    double c0 = calculate_c0(data.G0, data.rho, data.A0); // [cm/s]

    // in meter
    double C_e = (1e-2 * data.length) * (1e-4 * data.A0) / (1060 * std::pow(1e-2 * c0, 2));

    total_C_edge += C_e;
  }

  MPI_Allreduce(&total_C_edge, &total_C_edge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_C_edge;
}

double get_total_edge_capacitance(const std::vector<std::shared_ptr<GraphStorage>> &list) {
  double sum = 0;
  for (auto graph : list)
    sum += get_total_edge_capacitance(graph);
  return sum;
}

double get_total_flow(const std::vector<FlowData> &flows) {
  double sum = 0;
  for (auto &d : flows)
    sum += d.total_flow;
  return sum;
}

RCREstimator::RCREstimator(std::vector<std::shared_ptr<GraphStorage>> graph_list)
    : d_graph_list(std::move(graph_list)),
      d_total_C_edge(0),
      d_total_R(1.34),   // 1.34e8 Pa s / m^3   ([Pa s / m^-3] = 10^{-8} [kg/s cm])
      d_total_C(9.45e-1) // 9.e5e-9 m^3/Pa ([m^3 / Pa] = 10^8 [cm s^2 / kg])
{
  reset();
}

void RCREstimator::reset() {
  d_total_C_edge = get_total_edge_capacitance(d_graph_list);
}

void RCREstimator::set_total_C(double C) { d_total_C = C; }

void RCREstimator::set_total_R(double R) { d_total_R = R; }

std::map<size_t, RCRData> RCREstimator::estimate_parameters(const std::map<size_t, double> &flows, double total_flow) {
  std::map<size_t, RCRData> resistances;

  for (auto it : flows) {
    auto z = it.second / total_flow;
    double local_R = d_total_R / z;
    double C = z * (d_total_C - d_total_C_edge);

    resistances[it.first] = {local_R, C};
  }

  return resistances;
};

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
