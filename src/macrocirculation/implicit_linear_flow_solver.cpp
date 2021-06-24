////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "implicit_linear_flow_solver.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "interpolate_to_vertices.hpp"
#include "petsc.h"
#include "petsc/petsc_ksp.hpp"
#include "petsc/petsc_mat.hpp"
#include "petsc/petsc_vec.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

void assemble_mass(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, PetscVec &mass_vec) {
  // the squared norm of the legendre polynomials on [-1, +1]
  std::vector<double> legendre_weight = {2, 2. / 3., 0.4, 0.285714};

  std::vector<std::size_t> dof_indices(4, 0);

  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto edge = graph.get_edge(e_id);
    const auto local_dof_map = dof_map.get_local_dof_map(*edge);

    const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

    // assert that we have precalculated enough weights
    assert(num_basis_functions <= legendre_weight.size());

    dof_indices.resize(num_basis_functions);

    const auto edge_length = edge->get_physical_data().length;
    const auto micro_edge_length = edge_length / local_dof_map.num_micro_edges();
    const double edge_weight = micro_edge_length / 2;

    for (std::size_t component = 0; component < local_dof_map.num_components(); component += 1) {
      for (std::size_t local_micro_edge_id = 0; local_micro_edge_id < local_dof_map.num_micro_edges();
           local_micro_edge_id += 1) {
        local_dof_map.dof_indices(local_micro_edge_id, component, dof_indices);

        for (std::size_t i = 0; i < num_basis_functions; i += 1) {
          const double mass = legendre_weight[i] * edge_weight;
          mass_vec.set(static_cast<PetscInt>(dof_indices[i]), mass);
        }
      }
    }
  }

  for (const auto &v_id : graph.get_active_vertex_ids(mpi::rank(comm))) {
    auto &vertex = *graph.get_vertex(v_id);

    if (vertex.is_windkessel_outflow() || vertex.is_vessel_tree_outflow()) {
      auto &vertex_dof_map = dof_map.get_local_dof_map(vertex);
      auto &indices = vertex_dof_map.dof_indices();
      for (auto i : indices)
        mass_vec.set(static_cast<PetscInt>(i), 1.);
    }
  }
}

// TODO: Move somewhere else!!!
// TODO: Clean up from implicit advection solver
void interpolate_to_vertices(const MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMap &map,
                             const std::size_t component,
                             const PetscVec &dof_vector,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated) {
  points.clear();
  interpolated.clear();

  std::vector<std::size_t> dof_indices;
  std::vector<double> dof_vector_local;
  std::vector<double> evaluated_at_qps;

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    // we only write out embedded vessel segments
    if (!edge->has_embedding_data())
      continue;
    const auto &embedding = edge->get_embedding_data();

    auto local_dof_map = map.get_local_dof_map(*edge);

    if (embedding.points.size() == 2 && local_dof_map.num_micro_edges() > 1)
      linear_interpolate_points(embedding.points[0], embedding.points[1], local_dof_map.num_micro_edges(), points);
    else if (embedding.points.size() == local_dof_map.num_micro_edges() + 1)
      add_discontinuous_points(embedding.points, points);
    else
      throw std::runtime_error("this type of embedding is not implemented");

    FETypeNetwork fe(create_trapezoidal_rule(), local_dof_map.num_basis_functions() - 1);

    dof_indices.resize(local_dof_map.num_basis_functions());
    dof_vector_local.resize(local_dof_map.num_basis_functions());
    evaluated_at_qps.resize(fe.num_quad_points());

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    fe.reinit(h);

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
      extract_dof(dof_indices, dof_vector, dof_vector_local);
      const auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_vector_local);

      interpolated.push_back(boundary_values.left);
      interpolated.push_back(boundary_values.right);
    }
  }
}

LinearFlowSolver::LinearFlowSolver(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t degree)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      degree(degree),
      d_tau(0),
      u(std::make_shared<PetscVec>("u", *d_dof_map)),
      rhs(std::make_shared<PetscVec>("rhs", *d_dof_map)),
      A(std::make_shared<PetscMat>("A", *d_dof_map)),
      mass(std::make_shared<PetscVec>("mass", *d_dof_map)),
      linear_solver(std::make_shared<PetscKsp>(*A)) {
  assemble_mass(d_comm, *d_graph, *d_dof_map, *mass);
  u->zero();
  rhs->zero();
}

void LinearFlowSolver::setup(double tau) {
  d_tau = tau;
  assemble_matrix(tau);
}

void LinearFlowSolver::solve(double tau, double t) {
  if (std::abs(tau - d_tau) > 1e-14)
    throw std::runtime_error("changing time step size not supported");
  assemble_rhs(tau, t);
  linear_solver->solve(*rhs, *u);
}

Eigen::MatrixXd LinearFlowSolver::create_mass(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
  const auto &phi = fe.get_phi();
  const auto &JxW = fe.get_JxW();

  // TODO: This matrix is diagonal -> directly assemble it
  Eigen::MatrixXd m_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
    for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
      m_loc(j, i) = 0;
      for (int qp = 0; qp < phi[i].size(); qp += 1)
        m_loc(j, i) += phi[i][qp] * phi[j][qp] * JxW[qp];
    }
  }
  return m_loc;
}

Eigen::MatrixXd LinearFlowSolver::create_phi_grad_psi(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  Eigen::MatrixXd k_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
    for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
      k_loc(j, i) = 0;
      for (int qp = 0; qp < phi[i].size(); qp += 1)
        k_loc(j, i) += phi[i][qp] * dphi[j][qp] * JxW[qp];
    }
  }
  return k_loc;
}

enum class BoundaryPointType { Left,
                               Right };

Eigen::MatrixXd LinearFlowSolver::create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row, BoundaryPointType col) {
  Eigen::MatrixXd u_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  const auto left = [](size_t i) -> double { return std::pow(-1., i); };
  const auto right = [](size_t i) -> double { return 1.; };
  const auto phi = (col == BoundaryPointType::Left) ? left : right;
  const auto psi = (row == BoundaryPointType::Left) ? left : right;
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
    for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
      u_loc(j, i) = psi(j) * phi(i);
    }
  }
  return u_loc;
}

Eigen::MatrixXd LinearFlowSolver::create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType type) {
  Eigen::VectorXd u_loc(local_dof_map.num_basis_functions());
  const auto left = [](size_t i) -> double { return std::pow(-1., i); };
  const auto right = [](size_t i) -> double { return 1.; };
  const auto phi = (type == BoundaryPointType::Left) ? left : right;
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1)
    u_loc(j) = phi(j);
  return u_loc;
}

double LinearFlowSolver::get_C(const Edge &e) {
  const auto &data = e.get_physical_data();
  return data.A0 / (data.rho * std::pow(data.get_c0(), 2));
}

double LinearFlowSolver::get_L(const Edge &e) {
  const auto &data = e.get_physical_data();
  return data.rho / data.A0;
}

double LinearFlowSolver::get_R(const Edge &e) {
  const auto &data = e.get_physical_data();
  return 2 * (data.gamma + 2) * M_PI * data.viscosity / data.A0;
}

void LinearFlowSolver::assemble_matrix_cells(double tau) {
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);
    const auto &param = macro_edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    const double C = get_C(*macro_edge);
    const double L = get_L(*macro_edge);
    const double R = get_R(*macro_edge);

    std::vector<std::size_t> dof_indices_p(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_q(local_dof_map.num_basis_functions());

    const auto qf = create_gauss4();
    FETypeNetwork fe(qf, local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    auto k_loc{create_phi_grad_psi(fe, local_dof_map)};
    Eigen::MatrixXd k_pq = (-tau / C) * k_loc;
    Eigen::MatrixXd k_qp = (-tau / L) * k_loc;

    auto m_loc{create_mass(fe, local_dof_map)};
    Eigen::MatrixXd m_qq = tau * R * m_loc;

    for (const auto &edge : macro_edge->micro_edges()) {
      local_dof_map.dof_indices(edge, p_component, dof_indices_p);
      local_dof_map.dof_indices(edge, q_component, dof_indices_q);

      A->add(dof_indices_p, dof_indices_p, m_loc);
      A->add(dof_indices_p, dof_indices_q, k_pq);
      A->add(dof_indices_q, dof_indices_q, m_loc);
      A->add(dof_indices_q, dof_indices_p, k_qp);
      // A->add(dof_indices_q, dof_indices_q, m_qq);
    }
  }
}

void LinearFlowSolver::assemble_matrix_inner_boundaries(double tau) {
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);

    const double C = get_C(*macro_edge);
    const double L = get_L(*macro_edge);

    std::vector<std::size_t> dof_indices_p_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_q_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_p_right(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_q_right(local_dof_map.num_basis_functions());

    using BPT = BoundaryPointType;
    auto pattern_ll = create_boundary(local_dof_map, BPT::Left, BPT::Left);
    auto pattern_lr = create_boundary(local_dof_map, BPT::Left, BPT::Right);
    auto pattern_rl = create_boundary(local_dof_map, BPT::Right, BPT::Left);
    auto pattern_rr = create_boundary(local_dof_map, BPT::Right, BPT::Right);
    // + q^{up} phi(right)
    Eigen::MatrixXd u_pp_lr = tau * (1. / C) * (-0.5 * std::sqrt(C / L)) * pattern_rl;
    Eigen::MatrixXd u_pq_lr = tau * (1. / C) * (+0.5) * pattern_rl;
    Eigen::MatrixXd u_pp_ll = tau * (1. / C) * (+0.5 * std::sqrt(C / L)) * pattern_rr;
    Eigen::MatrixXd u_pq_ll = tau * (1. / C) * (+0.5) * pattern_rr;
    // - q^{up} phi(left)
    Eigen::MatrixXd u_pp_rr = tau * (1. / C) * (+0.5 * std::sqrt(C / L)) * pattern_ll;
    Eigen::MatrixXd u_pq_rr = tau * (1. / C) * (-0.5) * pattern_ll;
    Eigen::MatrixXd u_pp_rl = tau * (1. / C) * (-0.5 * std::sqrt(C / L)) * pattern_lr;
    Eigen::MatrixXd u_pq_rl = tau * (1. / C) * (-0.5) * pattern_lr;
    // - p^{up} phi(right)
    Eigen::MatrixXd u_qp_ll = tau * (1. / L) * (+0.5) * pattern_rr;
    Eigen::MatrixXd u_qq_ll = tau * (1. / L) * (+0.5 * std::sqrt(L / C)) * pattern_rr;
    Eigen::MatrixXd u_qp_lr = tau * (1. / L) * (+0.5) * pattern_rl;
    Eigen::MatrixXd u_qq_lr = tau * (1. / L) * (-0.5 * std::sqrt(L / C)) * pattern_rl;
    // + p^{up} phi(left)
    Eigen::MatrixXd u_qp_rl = tau * (1. / L) * (-0.5) * pattern_lr;
    Eigen::MatrixXd u_qq_rl = tau * (1. / L) * (-0.5 * std::sqrt(L / C)) * pattern_lr;
    Eigen::MatrixXd u_qp_rr = tau * (1. / L) * (-0.5) * pattern_ll;
    Eigen::MatrixXd u_qq_rr = tau * (1. / L) * (+0.5 * std::sqrt(L / C)) * pattern_ll;

    for (size_t micro_vertex_id = 1; micro_vertex_id < macro_edge->num_micro_vertices() - 1; micro_vertex_id += 1) {
      auto left_edge_id = micro_vertex_id - 1;
      auto right_edge_id = micro_vertex_id;

      local_dof_map.dof_indices(left_edge_id, p_component, dof_indices_p_left);
      local_dof_map.dof_indices(left_edge_id, q_component, dof_indices_q_left);
      local_dof_map.dof_indices(right_edge_id, p_component, dof_indices_p_right);
      local_dof_map.dof_indices(right_edge_id, q_component, dof_indices_q_right);

      // + q^{up} phi(right)
      A->add(dof_indices_p_left, dof_indices_p_right, u_pp_lr);
      A->add(dof_indices_p_left, dof_indices_q_right, u_pq_lr);
      A->add(dof_indices_p_left, dof_indices_p_left, u_pp_ll);
      A->add(dof_indices_p_left, dof_indices_q_left, u_pq_ll);
      // - q^{up} phi(left)
      A->add(dof_indices_p_right, dof_indices_p_right, u_pp_rr);
      A->add(dof_indices_p_right, dof_indices_q_right, u_pq_rr);
      A->add(dof_indices_p_right, dof_indices_p_left, u_pp_rl);
      A->add(dof_indices_p_right, dof_indices_q_left, u_pq_rl);
      // - p^{up} phi(right)
      A->add(dof_indices_q_left, dof_indices_p_left, u_qp_ll);
      A->add(dof_indices_q_left, dof_indices_q_left, u_qq_ll);
      A->add(dof_indices_q_left, dof_indices_p_right, u_qp_lr);
      A->add(dof_indices_q_left, dof_indices_q_right, u_qq_lr);
      // + p^{up} phi(left)
      A->add(dof_indices_q_right, dof_indices_p_left, u_qp_rl);
      A->add(dof_indices_q_right, dof_indices_q_left, u_qq_rl);
      A->add(dof_indices_q_right, dof_indices_p_right, u_qp_rr);
      A->add(dof_indices_q_right, dof_indices_q_right, u_qq_rr);
    }
  }
}

void LinearFlowSolver::assemble_rhs_inflow(double tau, double t) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);
    if (!vertex.is_inflow())
      continue;
    const auto q_in = vertex.get_inflow_value(t);
    auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
    const auto L = get_L(neighbor_edge);
    const auto C = get_C(neighbor_edge);
    const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;
    // b_p:
    std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
    std::vector<double> rhs_values_p(local_dof_map.num_basis_functions());
    for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
      // L^{-1} * tau * q_in(t) * phi_j(-1) = L^{-1} * tau * q_in(t) * (-1)^{j}
      rhs_values_p[j] = (-sigma / C) * tau * (-sigma * q_in) * std::pow(sigma, j);
    }
    rhs->add(dof_indices_p, rhs_values_p);
    // b_q:
    std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);
    std::vector<double> rhs_values_q(local_dof_map.num_basis_functions());
    for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
      // L^{-1} * tau * sqrt(L/C) * q_in(t) * phi_j(-1) = L^{-1} * tau * sqrt(L/C) * q_in(t) * (-1)^{j}
      rhs_values_q[j] = (-sigma / L) * tau * (sigma * std::sqrt(L / C)) * (sigma * q_in) * std::pow(sigma, j);
    }
    rhs->add(dof_indices_q, rhs_values_q);
  }
}

void LinearFlowSolver::assemble_matrix_inflow(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);
    if (!vertex.is_inflow())
      continue;
    auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
    std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
    std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
    const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;
    const auto L = get_L(neighbor_edge);
    const auto C = get_C(neighbor_edge);
    local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
    local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);
    auto pattern = neighbor_edge.is_pointing_to(v_idx) ? create_boundary(local_dof_map, BoundaryPointType::Right, BoundaryPointType::Right) : create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);
    Eigen::MatrixXd u_qp = sigma * tau * (1. / L) * pattern;
    Eigen::MatrixXd u_qq = tau * (1. / L) * (std::sqrt(L / C)) * pattern;
    A->add(dof_indices_q, dof_indices_p, u_qp);
    A->add(dof_indices_q, dof_indices_q, u_qq);
  }
}

void LinearFlowSolver::assemble_matrix_free_outflow(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);
    if (!vertex.is_free_outflow())
      continue;
    auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
    std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
    std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
    const auto L = get_L(neighbor_edge);
    const auto C = get_C(neighbor_edge);
    const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;
    local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
    local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);
    auto pattern = neighbor_edge.is_pointing_to(v_idx)
                     ? create_boundary(local_dof_map, BoundaryPointType::Right, BoundaryPointType::Right)
                     : create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);
    Eigen::MatrixXd u_pp = tau * (1. / C) * (0.5 * std::sqrt(C / L)) * pattern;
    Eigen::MatrixXd u_pq = tau * (1. / C) * sigma * (0.5) * pattern;
    Eigen::MatrixXd u_qp = tau * (1. / L) * sigma * (0.5) * pattern;
    Eigen::MatrixXd u_qq = tau * (1. / L) * (0.5 * std::sqrt(L / C)) * pattern;
    A->add(dof_indices_p, dof_indices_p, u_pp);
    A->add(dof_indices_p, dof_indices_q, u_pq);
    A->add(dof_indices_q, dof_indices_p, u_qp);
    A->add(dof_indices_q, dof_indices_q, u_qq);
  }
}

void LinearFlowSolver::assemble_matrix_nfurcations(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);
    if (!vertex.is_bifurcation())
      continue;
    const auto num_edges = vertex.get_edge_neighbors().size();

    // edge orientation: +1 if it point to the vertex, -1 else
    std::vector<double> sigma;
    std::vector<double> C;
    std::vector<double> L;
    std::vector<size_t> num_basis_functions;
    size_t all_basis_functions = 0;
    for (auto e_id : vertex.get_edge_neighbors()) {
      const auto &edge = *d_graph->get_edge(e_id);
      sigma.push_back(edge.is_pointing_to(v_idx) ? +1. : -1.);
      C.push_back(get_C(edge));
      L.push_back(get_L(edge));
      num_basis_functions.push_back(d_dof_map->get_local_dof_map(edge).num_basis_functions());
      all_basis_functions += num_basis_functions.back();
    }

    // inverse orientation matrix
    Eigen::MatrixXd orientation(2 * num_edges, 2 * num_edges);
    orientation.setZero();
    for (int k = 0; k < num_edges; k += 1) {
      orientation(2 * k, 2 * k) = tau * sigma[k] / L[k];
      orientation(2 * k + 1, 2 * k + 1) = tau * sigma[k] / C[k];
    }

    // contains the conditions on the upwind vectors
    Eigen::MatrixXd conditions(2 * num_edges, 2 * num_edges);
    conditions.setZero();
    // flow condition
    for (int k = 0; k < num_edges; k += 1)
      conditions(0, 2 * k + 1) = sigma[k];
    // pressure conditions
    for (int k = 1; k < num_edges; k += 1) {
      conditions(k, 0) = 1.;
      conditions(k, 2 * k) = -1.;
    }
    // characteristic conditions
    for (int k = 0; k < num_edges; k += 1) {
      conditions(num_edges + k, 2 * k) = sigma[k] * 0.5 * std::sqrt(C[k] / L[k]);
      conditions(num_edges + k, 2 * k + 1) = 0.5;
    }

    // inverse conditions
    Eigen::MatrixXd conditions_inverse = conditions.inverse();

    // evaluation matrix
    Eigen::MatrixXd evaluation_matrix(2 * num_edges, 2 * all_basis_functions);
    evaluation_matrix.setZero();
    size_t current_col = 0;
    for (size_t k = 0; k < num_edges; k += 1) {
      const auto e_id = vertex.get_edge_neighbors()[k];
      const auto &edge = *d_graph->get_edge(e_id);
      const double point = edge.is_pointing_to(v_idx) ? +1. : -1.;
      const auto psi = [=](size_t idx) { return std::pow(point, idx); };
      for (size_t i = 0; i < num_basis_functions[k]; i += 1) {
        evaluation_matrix(num_edges + k, current_col) = sigma[k] * 0.5 * std::sqrt(C[k] / L[k]) * psi(i);
        current_col += 1;
      }
      for (size_t i = 0; i < num_basis_functions[k]; i += 1) {
        evaluation_matrix(num_edges + k, current_col) = 0.5 * psi(i);
        current_col += 1;
      }
    }

    Eigen::MatrixXd mat = orientation * (conditions_inverse * evaluation_matrix);

    // apply test functions
    for (size_t j = 0; j < num_edges; j += 1) {
      const auto e_id_j = vertex.get_edge_neighbors()[j];
      const auto &edge_j = *d_graph->get_edge(e_id_j);

      // we only assemble rows for the given process.
      if (edge_j.rank() != mpi::rank(d_comm))
        continue;

      const double point = edge_j.is_pointing_to(v_idx) ? +1. : -1.;
      const auto psi = [=](size_t idx) { return std::pow(point, idx); };

      auto local_dof_map_j = d_dof_map->get_local_dof_map(edge_j);
      std::vector<size_t> dofs_p_j(local_dof_map_j.num_basis_functions());
      std::vector<size_t> dofs_q_j(local_dof_map_j.num_basis_functions());
      size_t micro_edge_idx_j = edge_j.is_pointing_to(v_idx) ? edge_j.num_micro_edges() - 1 : 0;
      local_dof_map_j.dof_indices(micro_edge_idx_j, p_component, dofs_p_j);
      local_dof_map_j.dof_indices(micro_edge_idx_j, q_component, dofs_q_j);

      for (size_t k = 0; k < num_edges; k += 1) {
        const auto e_id_k = vertex.get_edge_neighbors()[k];
        const auto &edge_k = *d_graph->get_edge(e_id_k);
        auto local_dof_map_k = d_dof_map->get_local_dof_map(edge_k);
        std::vector<size_t> dofs_p_k(local_dof_map_k.num_basis_functions());
        std::vector<size_t> dofs_q_k(local_dof_map_k.num_basis_functions());
        size_t micro_edge_idx_k = edge_k.is_pointing_to(v_idx) ? edge_k.num_micro_edges() - 1 : 0;
        local_dof_map_k.dof_indices(micro_edge_idx_k, p_component, dofs_p_k);
        local_dof_map_k.dof_indices(micro_edge_idx_k, q_component, dofs_q_k);

        Eigen::MatrixXd mat_el(local_dof_map_j.num_basis_functions(), local_dof_map_k.num_basis_functions());
        mat_el.setZero();

        // p^{up}:
        for (size_t row = 0; row < local_dof_map_j.num_basis_functions(); row += 1) {
          for (size_t col = 0; col < local_dof_map_j.num_basis_functions(); col += 1) {
            mat_el(row, col) = mat(2 * j, 2 * local_dof_map_k.num_basis_functions() * k + col) * psi(row);
          }
        }
        A->add(dofs_q_j, dofs_p_k, mat_el);
        for (size_t row = 0; row < local_dof_map_j.num_basis_functions(); row += 1) {
          for (size_t col = 0; col < local_dof_map_j.num_basis_functions(); col += 1) {
            mat_el(row, col) = mat(2 * j, 2 * local_dof_map_k.num_basis_functions() * k + local_dof_map_k.num_basis_functions() + col) * psi(row);
          }
        }
        A->add(dofs_q_j, dofs_q_k, mat_el);

        // q^{up}:
        for (size_t row = 0; row < local_dof_map_j.num_basis_functions(); row += 1) {
          for (size_t col = 0; col < local_dof_map_j.num_basis_functions(); col += 1) {
            mat_el(row, col) = mat(2 * j + 1, 2 * local_dof_map_k.num_basis_functions() * k + col) * psi(row);
          }
        }
        A->add(dofs_p_j, dofs_p_k, mat_el);
        for (size_t row = 0; row < local_dof_map_j.num_basis_functions(); row += 1) {
          for (size_t col = 0; col < local_dof_map_j.num_basis_functions(); col += 1) {
            mat_el(row, col) = mat(2 * j + 1, 2 * local_dof_map_k.num_basis_functions() * k + local_dof_map_k.num_basis_functions() + col) * psi(row);
          }
        }
        A->add(dofs_p_j, dofs_q_k, mat_el);
      }
    }
  }
}

void LinearFlowSolver::assemble_matrix_0d_model(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);

    if (vertex.is_windkessel_outflow() || vertex.is_vessel_tree_outflow()) {
      auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &local_dof_map_edge = d_dof_map->get_local_dof_map(neighbor_edge);
      auto &local_dof_map_vertex = d_dof_map->get_local_dof_map(vertex);

      std::vector<size_t> dof_indices_p(local_dof_map_edge.num_basis_functions());
      std::vector<size_t> dof_indices_q(local_dof_map_edge.num_basis_functions());

      auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
      const auto L = get_L(neighbor_edge);
      const auto C = get_C(neighbor_edge);
      const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;

      local_dof_map_edge.dof_indices(micro_edge_idx, p_component, dof_indices_p);
      local_dof_map_edge.dof_indices(micro_edge_idx, q_component, dof_indices_q);

      const auto &dof_indices_ptilde = local_dof_map_vertex.dof_indices();

      auto E = neighbor_edge.is_pointing_to(v_idx)
                 ? create_boundary(local_dof_map_edge, BoundaryPointType::Right, BoundaryPointType::Right)
                 : create_boundary(local_dof_map_edge, BoundaryPointType::Left, BoundaryPointType::Left);
      auto e = neighbor_edge.is_pointing_to(v_idx)
                 ? create_boundary(local_dof_map_edge, BoundaryPointType::Right)
                 : create_boundary(local_dof_map_edge, BoundaryPointType::Left);

      const double R0 = calculate_R1(neighbor_edge.get_physical_data());
      const auto R1 = vertex.is_windkessel_outflow() ? vertex.get_peripheral_vessel_data().resistance - R0 : vertex.get_vessel_tree_data().resistances[0];
      const auto C_tilde = vertex.is_windkessel_outflow() ? vertex.get_peripheral_vessel_data().compliance : vertex.get_vessel_tree_data().capacitances[0];

      const double alpha = sigma / (std::sqrt(C / L) + 1. / R0);

      Eigen::MatrixXd u_qp = (+sigma * tau / L) * alpha * (sigma * std::sqrt(C / L)) * E;
      Eigen::MatrixXd u_qq = (+sigma * tau / L) * alpha * E;
      Eigen::MatrixXd u_q_ptilde = (+sigma * tau / L) * alpha * sigma / R0 * e;

      Eigen::MatrixXd u_pp = (+sigma * tau / C) * sigma * (1. / R0 * alpha * sigma * std::sqrt(C / L)) * E;
      Eigen::MatrixXd u_pq = (+sigma * tau / C) * sigma * (1. / R0 * alpha) * E;
      Eigen::MatrixXd u_p_ptilde = (+sigma * tau / C) * sigma * (1. / R0 * (alpha * sigma / R0 - 1)) * e;

      Eigen::MatrixXd u_p0tilde_p0tilde(1, 1);
      u_p0tilde_p0tilde << 1. + tau / R0 / C_tilde + tau / R1 / C_tilde - tau / (R0 * C_tilde) * alpha * sigma / R0;
      Eigen::MatrixXd u_p0tilde_p = -tau / (R0 * C_tilde) * alpha * sigma * std::sqrt(C / L) * e.transpose();
      Eigen::MatrixXd u_p0tilde_q = -tau / (R0 * C_tilde) * alpha * e.transpose();

      A->add(dof_indices_q, dof_indices_q, u_qq);
      A->add(dof_indices_q, dof_indices_p, u_qp);
      A->add(dof_indices_q, {dof_indices_ptilde[0]}, u_q_ptilde);

      A->add(dof_indices_p, dof_indices_q, u_pq);
      A->add(dof_indices_p, dof_indices_p, u_pp);
      A->add(dof_indices_p, {dof_indices_ptilde[0]}, u_p_ptilde);

      // beginning
      A->add({dof_indices_ptilde[0]}, dof_indices_q, u_p0tilde_q);
      A->add({dof_indices_ptilde[0]}, dof_indices_p, u_p0tilde_p);
      A->add({dof_indices_ptilde[0]}, {dof_indices_ptilde[0]}, u_p0tilde_p0tilde);

      if (vertex.is_vessel_tree_outflow()) {
        const auto &data = vertex.get_vessel_tree_data();
        const auto &R = data.resistances;
        const auto &C_tilde2 = data.capacitances;

        for (size_t k = 1; k < R.size(); k += 1) {
          Eigen::MatrixXd mat_k_km1(1, 1);
          mat_k_km1 << -tau / (2 * R[k - 1]) / C_tilde2[k];
          Eigen::MatrixXd mat_k_k(1, 1);
          mat_k_k << 1. + tau / (2 * R[k - 1] * C_tilde2[k]) + tau / (R[k] * C_tilde2[k]);
          Eigen::MatrixXd mat_km1_k(1, 1);
          mat_km1_k << -tau / (R[k - 1] * C_tilde2[k - 1]);

          A->add({dof_indices_ptilde[k]}, {dof_indices_ptilde[k - 1]}, mat_k_km1);
          A->add({dof_indices_ptilde[k]}, {dof_indices_ptilde[k]}, mat_k_k);
          A->add({dof_indices_ptilde[k - 1]}, {dof_indices_ptilde[k]}, mat_km1_k);
        }
      }
    }
  }
}

void LinearFlowSolver::assemble_rhs_0d_model(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);

    if (vertex.is_windkessel_outflow() || vertex.is_vessel_tree_outflow()) {
      auto &edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);

      auto &local_dof_map_vertex = d_dof_map->get_local_dof_map(vertex);

      const auto &dof_indices_ptilde = local_dof_map_vertex.dof_indices();

      const double R0 = calculate_R1(edge.get_physical_data());
      const auto R1 = vertex.is_windkessel_outflow() ? vertex.get_peripheral_vessel_data().resistance - R0 : vertex.get_vessel_tree_data().resistances.back();
      const auto C_tilde = vertex.is_windkessel_outflow() ? vertex.get_peripheral_vessel_data().compliance : vertex.get_vessel_tree_data().capacitances.back();
      const auto p_out = vertex.is_windkessel_outflow() ? vertex.get_peripheral_vessel_data().p_out : vertex.get_vessel_tree_data().p_out;

      std::vector<double> value{tau * p_out / (R1 * C_tilde)};

      rhs->add({dof_indices_ptilde.back()}, value);
    }
  }
}

void LinearFlowSolver::assemble_matrix_characteristic(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);

    if (vertex.is_nonlinear_characteristic_inflow())
      throw std::runtime_error("linear solver found a nonlinear characteristic inflow boundary");

    if (vertex.is_linear_characteristic_inflow()) {
      auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);

      std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
      std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());

      auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
      const auto L_e = get_L(neighbor_edge);
      const auto C_e = get_C(neighbor_edge);
      const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;

      const auto &data = vertex.get_linear_characteristic_data();

      const auto L_v = data.L;
      const auto C_v = data.C;

      local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
      local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);

      auto E = neighbor_edge.is_pointing_to(v_idx)
                 ? create_boundary(local_dof_map, BoundaryPointType::Right, BoundaryPointType::Right)
                 : create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);

      const double beta_v = std::sqrt(C_v / L_v);
      const double beta_e = std::sqrt(C_e / L_e);
      const double alpha = 1. / (beta_v + beta_e);

      Eigen::MatrixXd u_pq = (-tau / C_e) * (1 - alpha * beta_e) * E;
      Eigen::MatrixXd u_pp = (-tau / C_e) * (-beta_e * (1 - alpha * beta_e)) * E;

      Eigen::MatrixXd u_qp = (-tau / L_e) * alpha * beta_e * E;
      Eigen::MatrixXd u_qq = (-tau / L_e) * alpha * (-1.) * E;

      A->add(dof_indices_p, dof_indices_q, u_pq);
      A->add(dof_indices_p, dof_indices_p, u_pp);
      A->add(dof_indices_q, dof_indices_q, u_qq);
      A->add(dof_indices_q, dof_indices_p, u_qp);
    }
  }
}

void LinearFlowSolver::assemble_rhs_characteristic(double tau) {
  for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_idx);

    if (vertex.is_nonlinear_characteristic_inflow())
      throw std::runtime_error("linear solver found a nonlinear characteristic inflow boundary");

    if (vertex.is_linear_characteristic_inflow()) {

      auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);

      std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
      std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());

      auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
      const auto L_e = get_L(neighbor_edge);
      const auto C_e = get_C(neighbor_edge);
      const double sigma = neighbor_edge.is_pointing_to(v_idx) ? +1 : -1;

      const auto &data = vertex.get_linear_characteristic_data();

      const auto L_v = data.L;
      const auto C_v = data.C;
      const auto q_v = data.q;
      const auto p_v = data.p;

      local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
      local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);

      const double beta_v = std::sqrt(C_v / L_v);
      const double beta_e = std::sqrt(C_e / L_e);
      const double alpha = 1. / (beta_v + beta_e);

      std::vector<double> rhs_values_p(local_dof_map.num_basis_functions());
      std::vector<double> rhs_values_q(local_dof_map.num_basis_functions());

      for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
        rhs_values_p[j] = (+tau/C_e) * (+beta_e*alpha*beta_v*p_v + beta_e *alpha*q_v) * std::pow(sigma, j);
        rhs_values_q[j] = (+tau/L_e) * alpha * ( beta_v * p_v + q_v) * std::pow(sigma, j);
      }

      rhs->add(dof_indices_p, rhs_values_p);
      rhs->add(dof_indices_q, rhs_values_q);
    }
  }
}

void LinearFlowSolver::assemble_matrix(double tau) {
  A->zero();
  assemble_matrix_cells(tau);
  assemble_matrix_inner_boundaries(tau);
  assemble_matrix_free_outflow(tau);
  assemble_matrix_inflow(tau);
  assemble_matrix_nfurcations(tau);
  assemble_matrix_0d_model(tau);
  assemble_matrix_characteristic(tau);
  A->assemble();
}

void LinearFlowSolver::assemble_rhs_cells() {
  CHKERRABORT(PETSC_COMM_WORLD, VecPointwiseMult(rhs->get_vec(), u->get_vec(), mass->get_vec()));
}

void LinearFlowSolver::assemble_rhs(double tau, double t) {
  rhs->zero();
  assemble_rhs_cells();
  assemble_rhs_inflow(tau, t);
  assemble_rhs_0d_model(tau);
  assemble_rhs_characteristic(tau);
  rhs->assemble();
}

void LinearFlowSolver::assemble(double tau, double t) {
  assemble_matrix(tau);
  assemble_rhs(tau, t);
}

} // namespace macrocirculation