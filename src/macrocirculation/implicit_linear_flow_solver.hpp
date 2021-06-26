////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_IMPLICIT_LINEAR_FLOW_SOLVER_HPP
#define TUMORMODELS_IMPLICIT_LINEAR_FLOW_SOLVER_HPP

#include "mpi.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class PetscVec;
class PetscMat;
class PetscKsp;
class FETypeNetwork;
class Point;
class LocalEdgeDofMap;
class Edge;
class Vertex;

void assemble_mass(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, PetscVec &mass_vec);

// TODO: Move somewhere else!!!
// TODO: Clean up from implicit advection solver
void interpolate_to_vertices(const MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMap &map,
                             const std::size_t component,
                             const PetscVec &dof_vector,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated);

class LinearFlowSolver {
public:
  LinearFlowSolver(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t degree);

  const size_t p_component = 0;
  const size_t q_component = 1;

  const PetscVec &get_solution() const { return *u; }

  void set_initial_value(double p, double q);

  void setup(double tau);

  void solve(double tau, double t);

  void assemble_matrix(double tau);

  void assemble_rhs(double tau, double t);

  void assemble(double tau, double t);

  static double get_C(const Edge &e);

  static double get_L(const Edge &e);

  static double get_R(const Edge &e);

  void get_1d_values_at_vertex(const Vertex& v, double& p, double& q) const;

private:
  static Eigen::MatrixXd create_mass(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map);

  static Eigen::MatrixXd create_phi_grad_psi(const FETypeNetwork &fe, const LocalEdgeDofMap &local_dof_map);

  enum class BoundaryPointType { Left,
                                 Right };

  static Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row, BoundaryPointType col);

  static Eigen::MatrixXd create_boundary(const LocalEdgeDofMap &local_dof_map, BoundaryPointType row);

  void assemble_matrix_cells(double tau);

  void assemble_matrix_inner_boundaries(double tau);

  void assemble_rhs_inflow(double tau, double t);

  void assemble_matrix_inflow(double tau);

  void assemble_matrix_free_outflow(double tau);

  void assemble_matrix_nfurcations(double tau);

  void assemble_rhs_cells();

  void assemble_matrix_0d_model(double tau);

  void assemble_rhs_0d_model(double tau);

  void assemble_matrix_characteristic(double tau);

  void assemble_rhs_characteristic(double tau);

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

#endif //TUMORMODELS_IMPLICIT_LINEAR_FLOW_SOLVER_HPP
