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
#include <macrocirculation/communication/mpi.hpp>
#include <macrocirculation/dof_map.hpp>
#include <macrocirculation/fe_type.hpp>
#include <memory>
#include <petsc.h>

#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_advection_solver.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/petsc/petsc_mat.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

namespace macrocirculation {

class LinearFlowSolver {
public:
  explicit LinearFlowSolver(std::shared_ptr<GraphStorage> graph, size_t degree)
      : d_graph(graph),
        d_dof_map(std::make_shared<DofMap>(d_graph->num_vertices(), d_graph->num_edges())),
        degree(degree),
        u(std::make_shared<PetscVec>("u", *d_dof_map)),
        rhs(std::make_shared<PetscVec>("rhs", *d_dof_map)),
        A(std::make_shared<PetscMat>("A", *d_dof_map)),
        linear_solver(std::make_shared<PetscKsp>(*A)) {
    d_dof_map->create(PETSC_COMM_WORLD, *d_graph, 2, degree, true);
  }

  const size_t p_component = 0;
  const size_t q_component = 1;

  void solve(double tau) {
  }

  void assemble() {
    A->zero();
    rhs->zero();

    for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(PETSC_COMM_WORLD))) {
      const auto macro_edge = d_graph->get_edge(e_id);

      const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);
      const auto &param = macro_edge->get_physical_data();
      const double h = param.length / local_dof_map.num_micro_edges();

      const double C = 1.;
      const double L = 1.;
      const double R = 1.;

      std::vector<std::size_t> dof_indices_p(local_dof_map.num_basis_functions());
      std::vector<std::size_t> dof_indices_q(local_dof_map.num_basis_functions());

      const auto qf = create_gauss4();
      FETypeNetwork fe(qf, local_dof_map.num_basis_functions() - 1);
      fe.reinit(h);

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

      auto k_pq = 1./C * k_loc;
      auto k_qp = 1./L * k_loc;

      // TODO: This matrix is diagonal -> directly assemble it
      Eigen::MatrixXd m_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
      for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
        for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
          m_loc(j, i) = 0;
          for (int qp = 0; qp < phi[i].size(); qp += 1)
            m_loc(j, i) += phi[i][qp] * phi[j][qp] * JxW[qp];
        }
      }

      auto m_qp = - R * m_loc;

      for (const auto &edge : macro_edge->micro_edges()) {
        local_dof_map.dof_indices(edge, p_component, dof_indices_p);
        local_dof_map.dof_indices(edge, q_component, dof_indices_q);


      }

      Eigen::MatrixXd u_pp(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
      Eigen::MatrixXd u_pq(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
      Eigen::MatrixXd u_qp(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
      Eigen::MatrixXd u_qq(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());


    }
  }


private:
  MPI_Comm d_comm{};

  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<DofMap> d_dof_map;

  size_t degree;

  std::shared_ptr<PetscVec> u;
  std::shared_ptr<PetscVec> rhs;
  std::shared_ptr<PetscMat> A;

  std::shared_ptr<PetscKsp> linear_solver;
};

} // namespace macrocirculation

int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 100;

  // initialize petsc
  MPI_Init(&argc, &argv);
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves advection problem"));

  std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

  const double velocity = 1;
  const double tau = 0.001;
  const double t_end = 2;
  const auto inflow_boundary_value = [](double t_now) -> double { return std::sin(M_PI * 3 * t_now); };
  // const auto inflow_boundary_value = [](double t_now) -> double { return 1; };

  // create the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  auto v0 = graph->create_vertex();
  auto v1 = graph->create_vertex();
  auto edge1 = graph->connect(*v0, *v1, num_micro_edges);

  edge1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1, 0, 0)}});
  edge1->add_physical_data({0, 0, 0, 1., 0, 0, 0.5});

  mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

  mc::LinearFlowSolver solver(graph, degree);

  solver.solve(tau);

  CHKERRQ(PetscFinalize());
}