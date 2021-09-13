////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "nonlinear_linear_coupling.hpp"

#include "communication/mpi.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

NonlinearLinearCoupling::
  NonlinearLinearCoupling(
    MPI_Comm comm,
    std::shared_ptr<GraphStorage> graph_nl,
    std::shared_ptr<GraphStorage> graph_li)
    : d_graph_nl(std::move(graph_nl)),
      d_graph_li(std::move(graph_li)),
      d_comm(comm),
      d_buffer_system(comm, 1) {}

void NonlinearLinearCoupling::add_coupled_vertices(const std::string &name) {
  add_coupled_vertices(name, name);
}

// TODO: calling this after having assembled the matrices will make the solvers fail
//       -> find a way to prevent this.
void NonlinearLinearCoupling::add_coupled_vertices(const std::string &name_nl, const std::string &name_li) {
  auto v_nl = d_graph_nl->find_vertex_by_name(name_nl);
  auto v_li = d_graph_li->find_vertex_by_name(name_li);
  assert(v_nl.is_leaf());
  assert(v_li.is_leaf());
  auto e_nl = d_graph_nl->get_edge(v_nl->get_edge_neighbors()[0]);
  auto e_li = d_graph_li->get_edge(v_li->get_edge_neighbors()[0]);
  const auto &data_nl = e_nl->get_physical_data();
  const auto &data_li = e_li->get_physical_data();
  v_nl->set_to_nonlinear_characteristic_inflow(data_li.G0, data_li.A0, data_li.rho, e_li->is_pointing_to(v_li->get_id()), 0, 0);
  v_li->set_to_linear_characteristic_inflow(linear::get_C(data_nl), linear::get_L(data_nl), e_nl->is_pointing_to(v_nl->get_id()), 0, 0);
  coupled_vertices.push_back({v_nl->get_id(), v_li->get_id()});
  Vertex::connect(d_graph_nl, *v_nl, d_graph_li, *v_li);
}

int NonlinearLinearCoupling::get_rank(GraphStorage &graph, Vertex &v) {
  assert(v.is_leaf());
  return graph.get_edge(v.get_edge_neighbors()[0])->rank();
}

void NonlinearLinearCoupling::update_linear_solver(const ExplicitNonlinearFlowSolver &nonlinear_solver, ImplicitLinearFlowSolver &linear_solver) {
  // send data
  for (auto vertex_pair : coupled_vertices) {
    auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
    auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

    auto sender = get_rank(*d_graph_nl, *v_nl);
    auto receiver = get_rank(*d_graph_li, *v_li);

    if (sender == mpi::rank(d_comm)) {
      double p, q;
      nonlinear_solver.get_1d_pq_values_at_vertex(*v_nl, p, q);
      d_buffer_system.get_send_buffer(receiver) << p << q;
      // we update our own characteristic inflow:
      v_li->update_linear_characteristic_inflow(p, q);
    }
  }

  d_buffer_system.start_communication();
  d_buffer_system.end_communication();

  for (auto vertex_pair : coupled_vertices) {
    auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
    auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

    auto sender = get_rank(*d_graph_nl, *v_nl);
    auto receiver = get_rank(*d_graph_li, *v_li);

    if (receiver == mpi::rank(d_comm)) {
      double p, q;
      d_buffer_system.get_receive_buffer(sender) >> p >> q;
      v_li->update_linear_characteristic_inflow(p, q);
    }
  }

  d_buffer_system.clear();
}

void NonlinearLinearCoupling::update_nonlinear_solver(const ImplicitLinearFlowSolver &linear_solver, ExplicitNonlinearFlowSolver &nonlinear_solver) {
  // send data
  for (auto vertex_pair : coupled_vertices) {
    auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
    auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

    auto sender = get_rank(*d_graph_li, *v_li);
    auto receiver = get_rank(*d_graph_nl, *v_nl);

    if (sender == mpi::rank(d_comm)) {
      double p, q;
      linear_solver.get_1d_pq_values_at_vertex(*v_li, p, q);
      d_buffer_system.get_send_buffer(receiver) << p << q;
      // we update our own nonlinear characteristic inflow:
      v_nl->update_nonlinear_characteristic_inflow(p, q);
    }
  }

  d_buffer_system.start_communication();
  d_buffer_system.end_communication();

  for (auto vertex_pair : coupled_vertices) {
    auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
    auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

    auto sender = get_rank(*d_graph_li, *v_li);
    auto receiver = get_rank(*d_graph_nl, *v_nl);

    if (receiver == mpi::rank(d_comm)) {
      double p, q;
      d_buffer_system.get_receive_buffer(sender) >> p >> q;
      v_nl->update_nonlinear_characteristic_inflow(p, q);
    }
  }

  d_buffer_system.clear();
}

} // namespace macrocirculation