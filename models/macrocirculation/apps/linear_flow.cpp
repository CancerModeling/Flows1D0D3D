////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <Eigen/Dense>
#include <cmath>
#include <gmm.h>
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_advection_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/petsc/petsc_mat.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

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

    if (!vertex.is_windkessel_outflow())
      continue;

    auto &vertex_dof_map = dof_map.get_local_dof_map(vertex);
    auto &indices = vertex_dof_map.dof_indices();
    for (auto i : indices)
      mass_vec.set(static_cast<PetscInt>(dof_indices[i]), 1.);
  }
}

// TODO: Move somewhere else!!!
// TODO: Clean up from implicit advection solver
void extract_dof(const std::vector<std::size_t> &dof_indices,
                 const PetscVec &global,
                 std::vector<double> &local) {
  assert(dof_indices.size() == local.size());

  for (std::size_t i = 0; i < dof_indices.size(); i += 1)
    local[i] = global.get(dof_indices[i]);
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

class LinearFlowSolver {
public:
  explicit LinearFlowSolver(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t degree)
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

  const size_t p_component = 0;
  const size_t q_component = 1;

  const PetscVec &get_solution() const { return *u; }

  void setup(double tau) {
    d_tau = tau;
    assemble_matrix(tau);
  }

  void solve(double tau, double t) {
    if (std::abs(tau - d_tau) > 1e-14)
      throw std::runtime_error("changing time step size not supported");
    assemble_rhs(tau, t);
    linear_solver->solve(*rhs, *u);
  }

  static Eigen::MatrixXd create_mass(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
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

  static Eigen::MatrixXd create_phi_grad_psi(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map) {
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

  static Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row, BoundaryPointType col) {
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

  static double get_C(const Edge &e) {
    const auto &data = e.get_physical_data();
    return data.A0 / (data.rho * std::pow(data.get_c0(), 2));
  }

  static double get_L(const Edge &e) {
    const auto &data = e.get_physical_data();
    return data.rho / data.A0;
  }

  static double get_R(const Edge &e) {
    const auto &data = e.get_physical_data();
    return 2 * (data.gamma + 2) * M_PI * data.viscosity / std::pow(data.A0, 2);
  }

  void assemble_matrix_cells(double tau) {
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
        A->add(dof_indices_q, dof_indices_q, m_qq);
      }
    }
  }

  void assemble_matrix_inner_boundaries(double tau) {
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

  void assemble_rhs_inflow(double tau, double t) {
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
      // b_p:
      std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
      local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
      std::vector<double> rhs_values_p(local_dof_map.num_basis_functions());
      for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
        // L^{-1} * tau * q_in(t) * phi_j(-1) = L^{-1} * tau * q_in(t) * (-1)^{j}
        rhs_values_p[j] = (1. / C) * tau * q_in * std::pow(-1, j);
      }
      rhs->add(dof_indices_p, rhs_values_p);
      // b_q:
      std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());
      local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);
      std::vector<double> rhs_values_q(local_dof_map.num_basis_functions());
      for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
        // L^{-1} * tau * sqrt(L/C) * q_in(t) * phi_j(-1) = L^{-1} * tau * sqrt(L/C) * q_in(t) * (-1)^{j}
        rhs_values_q[j] = (1. / L) * tau * std::sqrt(L / C) * q_in * std::pow(-1, j);
      }
      rhs->add(dof_indices_q, rhs_values_q);
    }
  }

  void assemble_matrix_inflow(double tau) {
    for (auto v_idx : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
      auto &vertex = *d_graph->get_vertex(v_idx);
      if (!vertex.is_inflow())
        continue;
      auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
      auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
      std::vector<size_t> dof_indices_p(local_dof_map.num_basis_functions());
      std::vector<size_t> dof_indices_q(local_dof_map.num_basis_functions());
      auto micro_edge_idx = neighbor_edge.is_pointing_to(v_idx) ? neighbor_edge.num_micro_edges() - 1 : 0;
      const auto L = get_L(neighbor_edge);
      const auto C = get_C(neighbor_edge);
      local_dof_map.dof_indices(micro_edge_idx, p_component, dof_indices_p);
      local_dof_map.dof_indices(micro_edge_idx, q_component, dof_indices_q);
      auto pattern_ll = create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);
      Eigen::MatrixXd u_qp = tau * (1. / L) * (-1.) * pattern_ll;
      Eigen::MatrixXd u_qq = tau * (1. / L) * (std::sqrt(L / C)) * pattern_ll;
      A->add(dof_indices_q, dof_indices_p, u_qp);
      A->add(dof_indices_q, dof_indices_q, u_qq);
    }
  }

  void assemble_matrix_free_outflow(double tau) {
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

  void assemble_matrix_nfurcations(double tau) {
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
        const double point = edge_j.is_pointing_to(v_idx) ? +1. : -1.;
        const auto psi = [=](size_t idx) { return std::pow(point, idx); };

        // TODO: Check for active edge!

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

  void assemble_matrix(double tau) {
    A->zero();
    assemble_matrix_cells(tau);
    assemble_matrix_inner_boundaries(tau);
    assemble_matrix_free_outflow(tau);
    assemble_matrix_inflow(tau);
    assemble_matrix_nfurcations(tau);
    A->assemble();
  }

  void assemble_rhs_cells() {
    CHKERRABORT(PETSC_COMM_WORLD, VecPointwiseMult(rhs->get_vec(), u->get_vec(), mass->get_vec()));
  }

  void assemble_rhs(double tau, double t) {
    rhs->zero();
    assemble_rhs_cells();
    assemble_rhs_inflow(tau, t);
    rhs->assemble();
  }

  void assemble(double tau, double t) {
    assemble_matrix(tau);
    assemble_rhs(tau, t);
  }

private:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<DofMap> d_dof_map;

  size_t degree;

  double d_tau;

  std::shared_ptr<PetscVec> u;
  std::shared_ptr<PetscVec> rhs;
  std::shared_ptr<PetscMat> A;

  std::shared_ptr<PetscVec> mass;

  std::shared_ptr<PetscKsp> linear_solver;
};

} // namespace macrocirculation

int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 50;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const double tau = 0.001;
    const double t_end = 0.5;

    const size_t output_interval = 1;

    // vessel parameters
    const double vessel_length = 42.2;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double gamma = 9;
    const double density = 1.028e-3;

    // create the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();

    auto v0 = graph->create_vertex();
    auto v1 = graph->create_vertex();
    auto v2 = graph->create_vertex();
    auto v3 = graph->create_vertex();
    auto v4 = graph->create_vertex();
    auto edge1 = graph->connect(*v0, *v1, num_micro_edges);
    auto edge2 = graph->connect(*v1, *v2, num_micro_edges);
    auto edge3 = graph->connect(*v1, *v3, num_micro_edges);
    // auto edge4 = graph->connect(*v1, *v4, num_micro_edges);
    auto edge4 = graph->connect(*v4, *v1, num_micro_edges);

    v0->set_to_inflow(mc::heart_beat_inflow(4.));
    v2->set_to_free_outflow();
    v3->set_to_free_outflow();
    v4->set_to_free_outflow();

    auto physical_data = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);

    edge1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1, 0, 0)}});
    edge1->add_physical_data(physical_data);
    edge2->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1, 1, 0)}});
    edge2->add_physical_data(physical_data);
    edge3->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1, -1, 0)}});
    edge3->add_physical_data(physical_data);
    // edge4->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(2, 0, 0)}});
    edge4->add_embedding_data({{ mc::Point(2, 0, 0), mc::Point(1, 0, 0) }});
    edge4->add_physical_data(physical_data);

    mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

    auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map->create(PETSC_COMM_WORLD, *graph, 2, degree, true);

    mc::LinearFlowSolver solver(PETSC_COMM_WORLD, graph, dof_map, degree);
    solver.setup(tau);

    mc::GraphPVDWriter writer(PETSC_COMM_WORLD, "./output", "linear_flow");

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver.solve(tau, t);

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;

        std::vector<mc::Point> points;
        std::vector<double> p_vertex_values;
        std::vector<double> q_vertex_values;
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.p_component, solver.get_solution(), points, p_vertex_values);
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.q_component, solver.get_solution(), points, q_vertex_values);

        writer.set_points(points);
        writer.add_vertex_data("p", p_vertex_values);
        writer.add_vertex_data("q", q_vertex_values);
        writer.write(t);
      }
    }
  }

  CHKERRQ(PetscFinalize());
}