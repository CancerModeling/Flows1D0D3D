////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TRANSPORT_H
#define TUMORMODELS_TRANSPORT_H

#include <memory>
#include <vector>
#include <cmath>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "gmm.h"

namespace macrocirculation {


class Transport {
public:
  Transport( MPI_Comm comm, std::shared_ptr< GraphStorage > graph, std::shared_ptr< DofMap > dof_map_flow, std::shared_ptr< DofMap > dof_map_transport )
    : d_comm(comm),
      d_graph(d_graph),
      d_dof_map_flow(d_dof_map_flow),
      d_dof_map_transport( d_dof_map_transport)
  {}

  void evaluate_macro_edge_boundary_values(const std::vector<double> &u_prev, const std::vector<double> &gamma_prev) {
    std::vector<std::size_t> dof_indices(4, 0);
    std::vector<double> local_dofs(4, 0);

    for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
      const auto edge = d_graph->get_edge(e_id);
      const auto &param = edge->get_physical_data();

      {
        const auto &local_dof_map = d_dof_map_flow->get_local_dof_map(*edge);
        const double h = param.length / double(local_dof_map.num_micro_edges());

        FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
        fe.reinit(h);

        dof_indices.resize(local_dof_map.num_basis_functions());
        local_dofs.resize(local_dof_map.num_basis_functions());

        local_dof_map.dof_indices(0, 0, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_Q_macro_edge_boundary_value_l[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

        local_dof_map.dof_indices(0, 1, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_A_macro_edge_boundary_value_l[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

        local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 0, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_Q_macro_edge_boundary_value_r[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).right;

        local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 1, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_A_macro_edge_boundary_value_r[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).right;
      }

      {
        const auto &local_dof_map = d_dof_map_transport->get_local_dof_map(*edge);
        const double h = param.length / double(local_dof_map.num_micro_edges());

        FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
        fe.reinit(h);

        dof_indices.resize(local_dof_map.num_basis_functions());
        local_dofs.resize(local_dof_map.num_basis_functions());

        local_dof_map.dof_indices(0, 0, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_gamma_macro_edge_boundary_value_l[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).left;

        local_dof_map.dof_indices(local_dof_map.num_micro_edges() - 1, 0, dof_indices);
        extract_dof(dof_indices, u_prev, local_dofs);
        d_gamma_macro_edge_boundary_value_r[edge->get_id()] = fe.evaluate_dof_at_boundary_points(local_dofs).right;
      }
    }
  }

  std::vector<double> &get_solution() { return d_solution; }


  void solve(double t, double dt, const std::vector< double >& u_old){
    evaluate_macro_edge_boundary_values( u_old, d_solution );
    // TODO: communication with mpi
    calculate_fluxes_at_nfurcations(t);
    assemble_rhs();
    apply_inverse_mass();
  }


private:

  void calculate_fluxes_at_nfurcations(double t)
  {
    for( auto& v_id : d_graph->get_vertex_ids() )
    {
      auto& vertex = *d_graph->get_vertex(v_id);

      if ( vertex.is_leaf() )
      {
        auto& edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);

        double Q = edge.is_pointing_to(vertex.get_id()) ? d_Q_macro_edge_boundary_value_r[edge.get_id()] : d_Q_macro_edge_boundary_value_l[edge.get_id()];
        double A = edge.is_pointing_to(vertex.get_id()) ? d_A_macro_edge_boundary_value_r[edge.get_id()] : d_A_macro_edge_boundary_value_l[edge.get_id()];

        double v = Q/A;

        const bool is_inflow = (v > 1e-8 && !edge.is_pointing_to(vertex.get_id())) || (v < 1e-8 && edge.is_pointing_to(vertex.get_id()));

        const double delta = 0.05;

        auto inflow_function = [=](double t) {
          if (vertex.is_inflow())
            return - 2*std::pow(t/delta, 3) + 3*std::pow(t/delta, 2);
          else
            return 0.;
        };

        if (is_inflow) {
          // inflow
          if (!edge.is_pointing_to(vertex.get_id())) {
            d_gamma_flux_l[edge.get_id()] = Q * inflow_function(t);
          }
          else
          {
            d_gamma_flux_r[edge.get_id()] = Q * inflow_function(t);
          }
        }
        else
        {
          // outflow:
          if (!edge.is_pointing_to(vertex.get_id())) {
            d_gamma_flux_l[edge.get_id()] = v*d_gamma_macro_edge_boundary_value_l[edge.get_id()];
          }
          else
          {
            d_gamma_flux_r[edge.get_id()] = v*d_gamma_macro_edge_boundary_value_r[edge.get_id()];
          }
        }
      }
      else if ( vertex.is_bifurcation() )
      {
        std::vector< double > Q;
        std::vector< double > A;
        std::vector< double > gamma;
        std::vector< bool > is_in;

        double Q_out = 0;
        double N_in = 0;

        for (auto e_id: vertex.get_edge_neighbors()) {
          auto& edge = *d_graph->get_edge(e_id);
          double Q_value = edge.is_pointing_to(vertex.get_id()) ? d_Q_macro_edge_boundary_value_r[edge.get_id()] : d_Q_macro_edge_boundary_value_l[edge.get_id()];
          double A_value = edge.is_pointing_to(vertex.get_id()) ? d_A_macro_edge_boundary_value_r[edge.get_id()] : d_A_macro_edge_boundary_value_l[edge.get_id()];
          double gamma_value = edge.is_pointing_to(vertex.get_id()) ? d_gamma_macro_edge_boundary_value_r[edge.get_id()] : d_gamma_macro_edge_boundary_value_l[edge.get_id()];
          double v = Q_value / A_value;
          const bool is_inflow_value = (v > 1e-8 && !edge.is_pointing_to(vertex.get_id())) || (v < 1e-8 && edge.is_pointing_to(vertex.get_id()));

          Q.push_back(Q_value);
          A.push_back(A_value);
          gamma.push_back(gamma_value);
          is_in.push_back((is_inflow_value));

          if (!is_inflow_value)
            Q_out += std::abs(Q_value);

          if (is_inflow_value)
            N_in += std::abs(Q_value) / A_value * gamma_value ;
        }

        // do nothing if there is no flow
        if ( gmm::vect_norminf(Q) < 1e-8 )
          continue;

        gmm::row_matrix<gmm::wsvector<double>> mat(vertex.get_edge_neighbors().size(), vertex.get_edge_neighbors().size());
        std::vector< double > rhs(vertex.get_edge_neighbors().size());

        for (std::size_t i = 0; i < vertex.get_edge_neighbors().size(); i+=1) {
          auto& edge = *d_graph->get_edge(vertex.get_edge_neighbors()[i]);
          double v = Q[edge.get_id()]/A[edge.get_id()];

          if (is_in[i]) {
            if (edge.is_pointing_to(v_id))
              d_gamma_flux_r[edge.get_id()] = v*d_gamma_macro_edge_boundary_value_r[edge.get_id()];
            else
              d_gamma_flux_l[edge.get_id()] = v*d_gamma_macro_edge_boundary_value_l[edge.get_id()];
          }
          else {
            double flux = v*N_in / Q_out * std::abs(Q[i])*A[i]/std::abs(Q[i]);
            if (edge.is_pointing_to(v_id))
              d_gamma_flux_r[edge.get_id()] = flux;
            else
              d_gamma_flux_l[edge.get_id()] = flux;
          }
        }
      }
      else
      {
        throw std::runtime_error("not implemented");
      }
    }
  }

  void apply_inverse_mass();

  MPI_Comm d_comm;

  std::vector< double > d_A_macro_edge_boundary_value_l;
  std::vector< double > d_A_macro_edge_boundary_value_r;
  std::vector< double > d_Q_macro_edge_boundary_value_l;
  std::vector< double > d_Q_macro_edge_boundary_value_r;
  std::vector< double > d_gamma_macro_edge_boundary_value_l;
  std::vector< double > d_gamma_macro_edge_boundary_value_r;

  std::vector< double > d_gamma_flux_l;
  std::vector< double > d_gamma_flux_r;

  std::vector< double > d_gamma_edge_l;
  std::vector< double > d_gamma_edge_r;
  std::vector< double > d_solution;

  std::shared_ptr< GraphStorage > d_graph;
  std::shared_ptr< DofMap > d_dof_map_flow;
  std::shared_ptr< DofMap > d_dof_map_transport;

};

}

#endif //TUMORMODELS_TRANSPORT_H
