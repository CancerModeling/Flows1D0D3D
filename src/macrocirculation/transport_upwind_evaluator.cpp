#include "transport_upwind_evaluator.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "implicit_transport_solver.hpp"
#include "explicit_transport_solver.hpp"
#include "vessel_formulas.hpp"

#include <cmath>
#include <utility>

namespace macrocirculation {

TransportUpwindEvaluator::TransportUpwindEvaluator(MPI_Comm comm,
                                                   std::shared_ptr<GraphStorage> graph,
                                                   std::shared_ptr<DofMap> dof_map,
                                                   std::shared_ptr<UpwindProvider> flow_upwind_provider)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_flow_upwind_provider(std::move(flow_upwind_provider)),
      d_gamma_boundary_evaluator(comm, d_graph, d_dof_map, 1),
      d_gamma_macro_edge_flux_l(d_graph->num_edges(), NAN),
      d_gamma_macro_edge_flux_r(d_graph->num_edges(), NAN),
      d_current_t(NAN) {}

void TransportUpwindEvaluator::init(double t, const PetscVec &u_prev) {
  d_current_t = t;

  d_gamma_boundary_evaluator.init(u_prev);

  calculate_nfurcation_fluxes(t, u_prev);
  calculate_inout_fluxes(t, u_prev);
}

void TransportUpwindEvaluator::set_inflow_function(std::function< double(double) > inflow_function){ d_inflow_function = std::move(inflow_function); }

void TransportUpwindEvaluator::calculate_nfurcation_fluxes(double t, const PetscVec &u_prev) {
  for (size_t vid : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    const auto &vertex = *d_graph->get_vertex(vid);

    std::vector<double> A_up(vertex.get_edge_neighbors().size());
    std::vector<double> Q_up(vertex.get_edge_neighbors().size());
    d_flow_upwind_provider->get_upwinded_values(t, vertex, A_up, Q_up);

    std::vector<double> sigmas = get_normals(*d_graph, vertex);

    std::vector<double> gammas(vertex.get_edge_neighbors().size());
    d_gamma_boundary_evaluator(vertex, gammas);

    if (vertex.is_bifurcation()) {
      std::vector<double> gamma_up(vertex.get_edge_neighbors().size());
      calculate_gamma_up_at_bifurcation(sigmas, Q_up, A_up, gammas, gamma_up);

      for (size_t k = 0; k < vertex.get_edge_neighbors().size(); k += 1) {
        auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[k]);
        if (edge.is_pointing_to(vid)) {
          d_gamma_macro_edge_flux_r[edge.get_id()] = gamma_up[k];
        } else {
          d_gamma_macro_edge_flux_l[edge.get_id()] = gamma_up[k];
        }
      }
    } else if (vertex.is_inflow_with_fixed_flow() || vertex.is_inflow_with_fixed_pressure()) {
      auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      const auto c_in = A_up[0] * d_inflow_function(t);
      if (edge.is_pointing_to(vid)) {
        d_gamma_macro_edge_flux_r[edge.get_id()] = c_in;
      } else {
        d_gamma_macro_edge_flux_l[edge.get_id()] = c_in;
      }
    // } else if (vertex.is_linear_characteristic_inflow()) {
    } else {
      auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      // const auto c_in = std::max(gammas.at(0), 1.);
      const auto c_in = gammas.at(0);
      if (edge.is_pointing_to(vid)) {
        d_gamma_macro_edge_flux_r[edge.get_id()] = c_in;
      } else {
        d_gamma_macro_edge_flux_l[edge.get_id()] = c_in;
      }
    }
  }
}

void TransportUpwindEvaluator::calculate_inout_fluxes(double t, const PetscVec &u_prev) {
  // TODO
}

bool TransportUpwindEvaluator::upwind_is_implemented(const Vertex& v) const
{
  return true;
  // return v.is_bifurcation() || v.is_inflow_with_fixed_pressure() || v.is_inflow_with_fixed_flow() || v.is_linear_characteristic_inflow();
}

void TransportUpwindEvaluator::get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &gamma_up) const {
  // evaluator was initialized with the correct time step
  if (d_current_t != t)
    throw std::runtime_error("FlowUpwindEvaluator was not initialized for the given time step");

  if (!upwind_is_implemented(v))
    throw std::runtime_error("not implemented");

  gamma_up.resize(v.get_edge_neighbors().size());

  for (size_t neighbor_edge_idx = 0; neighbor_edge_idx < v.get_edge_neighbors().size(); neighbor_edge_idx += 1) {
    const auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[neighbor_edge_idx]);

    if (edge.is_pointing_to(v.get_id())) {
      gamma_up[neighbor_edge_idx] = d_gamma_macro_edge_flux_r[edge.get_id()];
    } else {
      gamma_up[neighbor_edge_idx] = d_gamma_macro_edge_flux_l[edge.get_id()];
    }
  }
}

} // namespace macrocirculation