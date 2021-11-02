////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "mpi.h"
#include <cmath>

#include "linearized_flow_upwind_evaluator.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "upwinding_formulas_linearized_flow.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

template<typename VectorType>
inline void compute_inner_boundary_values_on_macro_edge(const Edge &edge,
                                                        const DofMap &dof_map,
                                                        const VectorType &u,
                                                        size_t component_index,
                                                        std::vector<double> &boundary_values_l,
                                                        std::vector<double> &boundary_values_r) {
  auto &local_dof_map = dof_map.get_local_dof_map(edge);

  assert(local_dof_map.num_micro_edges() == boundary_values_l.size());
  assert(local_dof_map.num_micro_edges() == boundary_values_r.size());

  const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

  const auto &param = edge.get_physical_data();
  const double h = param.length / local_dof_map.num_micro_edges();

  std::vector<std::size_t> dof_indices(num_basis_functions, 0);
  std::vector<double> dof_values(num_basis_functions, 0);

  FETypeNetwork fe(create_trapezoidal_rule(), num_basis_functions - 1);
  fe.reinit(h);

  for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
    local_dof_map.dof_indices(micro_edge_id, component_index, dof_indices);
    extract_dof(dof_indices, u, dof_values);
    auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_values);
    boundary_values_l[micro_edge_id] = boundary_values.left;
    boundary_values_r[micro_edge_id] = boundary_values.right;
  }
}

LinearizedFlowUpwindEvaluator::LinearizedFlowUpwindEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_q_boundary_evaluator(comm, d_graph, d_dof_map, 1),
      d_p_boundary_evaluator(comm, d_graph, d_dof_map, 0),
      d_q_macro_edge_flux_l(d_graph->num_edges(), NAN),
      d_q_macro_edge_flux_r(d_graph->num_edges(), NAN),
      d_p_macro_edge_flux_l(d_graph->num_edges(), NAN),
      d_p_macro_edge_flux_r(d_graph->num_edges(), NAN),
      d_current_t(NAN) {}

void LinearizedFlowUpwindEvaluator::init(double t, const std::vector<double> &u_prev) {
  init_generic(t, u_prev);
}

void LinearizedFlowUpwindEvaluator::init(double t, const PetscVec &u_prev) {
  init_generic(t, u_prev);
}

void LinearizedFlowUpwindEvaluator::get_fluxes_on_macro_edge(double t, const Edge &edge, const std::vector<double> &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
  get_fluxes_on_macro_edge_generic(t, edge, u_prev, p_up, q_up);
}

void LinearizedFlowUpwindEvaluator::get_fluxes_on_macro_edge(double t, const Edge &edge, const PetscVec &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
  get_fluxes_on_macro_edge_generic(t, edge, u_prev, p_up, q_up);
}

void LinearizedFlowUpwindEvaluator::get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &p_up, std::vector<double> &q_up) const {
  // evaluator was initialized with the correct time step
  if (d_current_t != t)
    throw std::runtime_error("LinearizedFlowUpwindEvaluator was not initialized for the given time step");

  p_up.resize(v.get_edge_neighbors().size());
  q_up.resize(v.get_edge_neighbors().size());

  for (size_t neighbor_edge_idx = 0; neighbor_edge_idx < v.get_edge_neighbors().size(); neighbor_edge_idx += 1) {
    const auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[neighbor_edge_idx]);

    if (edge.is_pointing_to(v.get_id())) {
      p_up[neighbor_edge_idx] = d_p_macro_edge_flux_r[edge.get_id()];
      q_up[neighbor_edge_idx] = d_q_macro_edge_flux_r[edge.get_id()];
    } else {
      p_up[neighbor_edge_idx] = d_p_macro_edge_flux_l[edge.get_id()];
      q_up[neighbor_edge_idx] = d_q_macro_edge_flux_l[edge.get_id()];
    }
  }

  //std::cout << t << " " << v.get_id() << " " << p_up << " " << q_up << std::endl;
  //std::cout << t << " " << v.get_id() << " q_l = " << d_q_macro_edge_flux_l << std::endl;
  //std::cout << t << " " << v.get_id() << " q_r = " << d_q_macro_edge_flux_r << std::endl;
}

template<typename VectorType>
void LinearizedFlowUpwindEvaluator::init_generic(double t, const VectorType &u_prev) {
  d_current_t = t;

  d_p_boundary_evaluator.init(u_prev);
  d_q_boundary_evaluator.init(u_prev);

  calculate_nfurcation_fluxes(u_prev);
  calculate_inout_fluxes(t, u_prev);
}

template<typename VectorType>
void LinearizedFlowUpwindEvaluator::get_fluxes_on_macro_edge_generic(double t, const Edge &edge, const VectorType &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
  // evaluator was initialized with the correct time step
  if (d_current_t != t)
    throw std::runtime_error("LinearizedFlowUpwindEvaluator was not initialized for the given time step");

  assert(p_up.size() == edge.num_micro_vertices());
  assert(q_up.size() == edge.num_micro_vertices());

  std::vector<double> p_l(edge.num_micro_edges(), 0);
  std::vector<double> p_r(edge.num_micro_edges(), 0);
  std::vector<double> q_l(edge.num_micro_edges(), 0);
  std::vector<double> q_r(edge.num_micro_edges(), 0);

  compute_inner_boundary_values_on_macro_edge(edge, *d_dof_map, u_prev, 0, p_l, p_r);
  compute_inner_boundary_values_on_macro_edge(edge, *d_dof_map, u_prev, 1, q_l, q_r);

  assert(edge.has_physical_data());
  auto param = edge.get_physical_data();

  // TODO: move this calculation into its own function
  const auto alpha = std::sqrt(linear::get_C(param) / linear::get_L(param));

  for (std::size_t micro_vertex_id = 1; micro_vertex_id < edge.num_micro_vertices() - 1; micro_vertex_id += 1) {
    linearized::inner_boundary(alpha, p_l[micro_vertex_id], q_l[micro_vertex_id], p_r[micro_vertex_id - 1], q_r[micro_vertex_id - 1], p_up[micro_vertex_id], q_up[micro_vertex_id]);
  }

  // update left fluxes
  p_up[0] = d_p_macro_edge_flux_l[edge.get_id()];
  q_up[0] = d_q_macro_edge_flux_l[edge.get_id()];

  // update right fluxes
  p_up[edge.num_micro_vertices() - 1] = d_p_macro_edge_flux_r[edge.get_id()];
  q_up[edge.num_micro_vertices() - 1] = d_q_macro_edge_flux_r[edge.get_id()];
}

template<typename VectorType>
void LinearizedFlowUpwindEvaluator::calculate_nfurcation_fluxes(const VectorType &u_prev) {
  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // we only handle bifurcations
    if (!vertex->is_bifurcation())
      continue;

    const size_t num_vessels = vertex->get_edge_neighbors().size();

    // get edges
    std::vector<std::shared_ptr<Edge>> e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      e.push_back(d_graph->get_edge(vertex->get_edge_neighbors()[vessel_idx]));

    // check orientation
    std::vector<double> sigma;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      sigma.push_back(e[vessel_idx]->is_pointing_to(vertex->get_id()) ? +1 : -1);

    // get data
    std::vector<VesselParameters> params;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      const auto &data_e = e[vessel_idx]->get_physical_data();
      params.emplace_back(data_e.G0, data_e.A0, data_e.rho);
    }

    std::vector<double> p_e(num_vessels);
    std::vector<double> q_e(num_vessels);
    d_p_boundary_evaluator(*vertex, p_e);
    d_q_boundary_evaluator(*vertex, q_e);

    std::vector<double> p_up(num_vessels, 0);
    std::vector<double> q_up(num_vessels, 0);

    linearized::nfurcation_boundary(p_e, q_e, params, sigma, p_up, q_up);

    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      if (e[vessel_idx]->is_pointing_to(v_id)) {
        d_p_macro_edge_flux_r[e[vessel_idx]->get_id()] = p_up[vessel_idx];
        d_q_macro_edge_flux_r[e[vessel_idx]->get_id()] = q_up[vessel_idx];
      } else {
        d_p_macro_edge_flux_l[e[vessel_idx]->get_id()] = p_up[vessel_idx];
        d_q_macro_edge_flux_l[e[vessel_idx]->get_id()] = q_up[vessel_idx];
      }
    }
  }
}

template<typename VectorType>
void LinearizedFlowUpwindEvaluator::calculate_inout_fluxes(double t, const VectorType &u_prev) {
  for (const auto &v_id : d_graph->get_active_and_connected_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // we are only interested in leaves and ignore the rest:
    if (!vertex->is_leaf())
      continue;

    const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
    const auto &param = edge->get_physical_data();
    const double sigma = edge->is_pointing_to(vertex->get_id()) ? +1. : -1.;

    std::vector<double> p_value;
    std::vector<double> q_value;
    d_p_boundary_evaluator(*vertex, p_value);
    d_q_boundary_evaluator(*vertex, q_value);

    // TODO: move this calculation into its own function
    const auto alpha = std::sqrt(linear::get_C(param) / linear::get_L(param));

    double p_up = NAN;
    double q_up = NAN;

    if (vertex->is_inflow_with_fixed_flow()) {
      q_up = -sigma * vertex->get_inflow_value(t);
      p_up = p_value[0] + 1. / alpha * (q_value[0] - q_up);
    } else if (vertex->is_inflow_with_fixed_pressure()) {
      p_up = vertex->get_inflow_value(t);
      q_up = sigma * alpha * p_value[0] + q_value[0] - sigma * alpha * p_up;
    } else if (vertex->is_free_outflow()) {
      p_up = 0.5 * (p_value[0] + sigma / alpha * q_value[0]);
      q_up = sigma * alpha * p_up;
    } else if (vertex->is_linear_characteristic_inflow()) {
      const auto C_tilde = vertex->get_linear_characteristic_data().C;
      const auto L_tilde = vertex->get_linear_characteristic_data().L;

      const auto p_tilde = vertex->get_linear_characteristic_data().p;
      const auto q_tilde = vertex->get_linear_characteristic_data().q;

      const auto C = linear::get_C(param);
      const auto L = linear::get_L(param);

      const auto beta = std::sqrt(C / L);
      const auto beta_tilde = std::sqrt(C_tilde / L_tilde);
      const auto gamma = beta + beta_tilde;

      p_up = 1. / gamma * (beta_tilde * p_tilde + beta * p_value[0] + q_tilde + sigma * q_value[0]);
      q_up = sigma * ( beta * p_value[0] + sigma * q_value[0] - beta * p_up);
    } else if (vertex->is_windkessel_outflow() || vertex->is_vessel_tree_outflow()) {
      const double R0 = calculate_R1(edge->get_physical_data());
      const double C = linear::get_C(edge->get_physical_data());
      const double L = linear::get_L(edge->get_physical_data());
      const double beta = std::sqrt(C/L);

      const auto& vdof_map = d_dof_map->get_local_dof_map(*vertex);
      const auto& vertex_dof_indices = vdof_map.dof_indices();

      std::vector< double > vertex_dof_values (vdof_map.dof_indices().size());
      extract_dof(vertex_dof_indices, u_prev, vertex_dof_values);

      const double p_c = vertex_dof_values[0];

      p_up = 1./(beta + 1./R0) * (beta * p_value[0] + sigma * q_value[0] + p_c / R0);
      q_up = sigma / R0 * (p_up - p_c);
    } else {
      throw std::runtime_error("boundary type not implemented yet in LinearizedFlowUpwindEvaluator.");
    }

    // write back into vectors:
    if (edge->is_pointing_to(vertex->get_id())) {
      d_q_macro_edge_flux_r[edge->get_id()] = q_up;
      d_p_macro_edge_flux_r[edge->get_id()] = p_up;
    } else {
      d_q_macro_edge_flux_l[edge->get_id()] = q_up;
      d_p_macro_edge_flux_l[edge->get_id()] = p_up;
    }
  }
}

} // namespace macrocirculation