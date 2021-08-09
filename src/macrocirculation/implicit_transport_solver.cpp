////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "implicit_transport_solver.hpp"

#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "petsc/petsc_ksp.hpp"
#include "petsc/petsc_vec.hpp"
#include "petsc_assembly_blocks.hpp"

namespace macrocirculation {

ConstantUpwindProvider::ConstantUpwindProvider(double speed)
    : d_speed(speed) {}

ConstantUpwindProvider::~ConstantUpwindProvider() = default;

void ConstantUpwindProvider::get_values_at_qp(double t,
                                              const Edge &edge,
                                              size_t micro_edge,
                                              const QuadratureFormula &qf,
                                              std::vector<double> &v_qp) const {
  assert(qf.size() == v_qp.size());

  std::uniform_real_distribution<double> distribution(0.5, 1.);

  for (size_t k = 0; k < qf.size(); k += 1)
    v_qp[k] = d_speed;
}

/*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
void ConstantUpwindProvider::get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const {
  assert(v_qp.size() == edge.num_micro_vertices());

  std::uniform_real_distribution<double> distribution(0.5, 1.);

  for (size_t k = 0; k < edge.num_micro_vertices(); k += 1)
    v_qp[k] = d_speed;
}

void ConstantUpwindProvider::get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const {
  assert(v.get_edge_neighbors().size() == 1);
  assert(A.size() == 1);
  assert(Q.size() == 1);

  std::uniform_real_distribution<double> distribution(0.5, 1.);

  A[0] = 1.;
  Q[0] = d_speed;
}

ImplicitTransportSolver::ImplicitTransportSolver(MPI_Comm comm,
                                                 std::shared_ptr<GraphStorage> graph,
                                                 std::shared_ptr<DofMap> dof_map,
                                                 std::shared_ptr<UpwindProvider> upwind_provider,
                                                 size_t degree)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_upwind_provider(std::move(upwind_provider)),
      d_degree(degree),
      u(std::make_shared<PetscVec>("u", *d_dof_map)),
      rhs(std::make_shared<PetscVec>("rhs", *d_dof_map)),
      A(std::make_shared<PetscMat>("A", *d_dof_map)),
      mass(std::make_shared<PetscVec>("mass", *d_dof_map)),
      linear_solver(PetscKsp::create_with_pc_ilu(*A)) {
  assemble_mass(d_comm, *d_graph, *d_dof_map, *mass);
  u->zero();
  rhs->zero();
  // initialize zero pattern of matrix with velocity 1 vector:
  // otherwise in some vessels the velocity might be zero and thus the sparsity pattern would change
  sparsity_pattern();
}

void ImplicitTransportSolver::solve(double tau, double t) {
  assemble(tau, t);
  linear_solver->solve(*rhs, *u);
}

void ImplicitTransportSolver::assemble_rhs_cells(double tau, double t, const UpwindProvider &upwind_provider) {
  CHKERRABORT(PETSC_COMM_WORLD, VecPointwiseMult(rhs->get_vec(), u->get_vec(), mass->get_vec()));
}

void ImplicitTransportSolver::assemble(double tau, double t) {
  assemble_matrix(tau, t, *d_upwind_provider);
  assemble_rhs(tau, t, *d_upwind_provider);
}

void ImplicitTransportSolver::sparsity_pattern()
{
  ConstantUpwindProvider upwind_provider_plus(1.);
  ConstantUpwindProvider upwind_provider_minus(-0.5);
  const double tau = 1e-3;
  const double t = 1.;
  A->zero();
  assemble_matrix_cells(tau, t, upwind_provider_plus);
  assemble_matrix_cells(tau, t, upwind_provider_minus);
  assemble_matrix_inner_boundaries(tau, t, upwind_provider_plus);
  assemble_matrix_inner_boundaries(tau, t, upwind_provider_minus);
  assemble_matrix_outflow(tau, t, upwind_provider_plus);
  assemble_matrix_outflow(tau, t, upwind_provider_minus);
  A->assemble();
}

void ImplicitTransportSolver::assemble_matrix(double tau, double t, const UpwindProvider &upwind_provider) {
  A->zero();
  assemble_matrix_cells(tau, t, upwind_provider);
  assemble_matrix_inner_boundaries(tau, t, upwind_provider);
  assemble_matrix_outflow(tau, t, upwind_provider);
  A->assemble();
}

void ImplicitTransportSolver::assemble_rhs(double tau, double t, const UpwindProvider &upwind_provider) {
  rhs->zero();
  assemble_rhs_cells(tau, t, upwind_provider);
  assemble_rhs_inflow(tau, t, upwind_provider);
  rhs->assemble();
}

Eigen::MatrixXd create_QA_phi_grad_psi(const FETypeNetwork &fe,
                                       const LocalEdgeDofMap &local_dof_map,
                                       const std::vector<double> &v_qp) {
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  Eigen::MatrixXd k_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
    for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
      k_loc(j, i) = 0;
      for (int qp = 0; qp < phi[i].size(); qp += 1)
        k_loc(j, i) += v_qp[qp] * phi[i][qp] * dphi[j][qp] * JxW[qp];
    }
  }
  return k_loc;
}

void ImplicitTransportSolver::assemble_matrix_cells(double tau, double t, const UpwindProvider &upwind_provider) {
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);
    const auto &param = macro_edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    std::vector<std::size_t> dof_indices_gamma(local_dof_map.num_basis_functions());

    const auto qf = create_gauss4();
    FETypeNetwork fe(qf, local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    std::vector<double> v_qp(qf.size());

    auto m_loc{create_mass(fe, local_dof_map)};

    for (const auto &edge : macro_edge->micro_edges()) {
      local_dof_map.dof_indices(edge, 0, dof_indices_gamma);

      upwind_provider.get_values_at_qp(t, *macro_edge, edge.get_local_id(), fe.get_quadrature_formula(), v_qp);

      auto k_loc{create_QA_phi_grad_psi(fe, local_dof_map, v_qp)};

      Eigen::MatrixXd mat = m_loc - tau * k_loc;
      A->add(dof_indices_gamma, dof_indices_gamma, mat);
    }
  }
}

void ImplicitTransportSolver::assemble_matrix_inner_boundaries(double tau, double t, const UpwindProvider &upwind_provider) {
  std::vector<double> v_up;

  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    v_up.resize(macro_edge->num_micro_vertices());

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);

    std::vector<std::size_t> dof_indices_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_right(local_dof_map.num_basis_functions());

    using BPT = BoundaryPointType;
    auto pattern_ll = create_boundary(local_dof_map, BPT::Left, BPT::Left);
    auto pattern_lr = create_boundary(local_dof_map, BPT::Left, BPT::Right);
    auto pattern_rl = create_boundary(local_dof_map, BPT::Right, BPT::Left);
    auto pattern_rr = create_boundary(local_dof_map, BPT::Right, BPT::Right);

    for (size_t micro_vertex_id = 1; micro_vertex_id < macro_edge->num_micro_vertices() - 1; micro_vertex_id += 1) {
      auto left_edge_id = micro_vertex_id - 1;
      auto right_edge_id = micro_vertex_id;

      local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);
      local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);

      upwind_provider.get_upwinded_values(t, *macro_edge, v_up);

      const double v = v_up[micro_vertex_id];

      if (v > 0) {
        Eigen::MatrixXd mat_ll = +tau * v * pattern_rr;
        Eigen::MatrixXd mat_rl = -tau * v * pattern_lr;

        A->add(dof_indices_left, dof_indices_left, mat_ll);
        A->add(dof_indices_right, dof_indices_left, mat_rl);
      } else {
        Eigen::MatrixXd mat_lr = +tau * v * pattern_rl;
        Eigen::MatrixXd mat_rr = -tau * v * pattern_ll;

        A->add(dof_indices_left, dof_indices_right, mat_lr);
        A->add(dof_indices_right, dof_indices_right, mat_rr);
      }
    }
  }
}

void ImplicitTransportSolver::assemble_rhs_inflow(double tau, double t, const UpwindProvider &upwind_provider) {
  std::vector<double> Q_up(1);
  std::vector<double> A_up(1);

  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_id);
    if (!vertex.is_inflow())
      continue;

    // TODO: generalize this to allow different functions
    const auto c_in = inflow_function(t);

    auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_id) ? neighbor_edge.num_micro_edges() - 1 : 0;
    const double sigma = neighbor_edge.is_pointing_to(v_id) ? +1 : -1;

    upwind_provider.get_upwinded_values(t, vertex, A_up, Q_up);

    const double v_up = Q_up[0] / A_up[0];

    std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, 0, dof_indices);
    std::vector<double> rhs_values(local_dof_map.num_basis_functions());
    for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
      // L^{-1} * tau * q_in(t) * phi_j(-1) = L^{-1} * tau * q_in(t) * (-1)^{j}
      rhs_values[j] = tau * (-sigma * v_up * c_in) * std::pow(sigma, j);
    }

    rhs->add(dof_indices, rhs_values);
  }
}

void ImplicitTransportSolver::assemble_matrix_outflow(double tau, double t, const UpwindProvider &upwind_provider) {
  std::vector<double> Q_up(1);
  std::vector<double> A_up(1);

  for (const auto &v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_id);

    if (!vertex.is_leaf())
      continue;
    if (vertex.is_inflow())
      continue;

    auto &neighbor_edge = *d_graph->get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = d_dof_map->get_local_dof_map(neighbor_edge);
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_id) ? neighbor_edge.num_micro_edges() - 1 : 0;

    std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, 0, dof_indices);

    const double sigma = neighbor_edge.is_pointing_to(v_id) ? +1 : -1;

    auto pattern = neighbor_edge.is_pointing_to(v_id)
                     ? create_boundary(local_dof_map, BoundaryPointType::Right, BoundaryPointType::Right)
                     : create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);

    upwind_provider.get_upwinded_values(t, vertex, A_up, Q_up);

    const double v_up = Q_up[0] / A_up[0];

    Eigen::MatrixXd mat = sigma * tau * v_up * pattern;

    A->add(dof_indices, dof_indices, mat);
  }
}

} // namespace macrocirculation
