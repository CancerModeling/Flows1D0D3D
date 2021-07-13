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

/*! @brief An implicit euler for the  linear flow equations. */
class ImplicitLinearFlowSolver {
public:
  ImplicitLinearFlowSolver(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, size_t degree);

  /*! @brief The component index of the pressure p inside the dof map.
   *         This can be used inside a local dof map to get the dof indices for a single component. */
  static const size_t p_component = 0;

  /*! @brief The component index of the flow q inside the dof map.
   *         This can be used inside a local dof map to get the dof indices for a single component. */
  static const size_t q_component = 1;

  const DofMap &get_dof_map() const;

  const PetscVec &get_solution() const { return *u; }

  /*! @brief Requests the Jacobi preconditioner for the linear solver, which yields (more) reproducible results for unit-testing. */
  void use_pc_jacobi();

  /*! @brief Sets the initial constant value of p and q. */
  void set_initial_value(double p, double q);

  /*! @brief Sets up the matrices for the given time step size tau. */
  void setup(double tau);

  /*! @brief Solves the system. For the _next_ time step t.
   *
   * @param tau The current time step width.
   * @param t The (future) time for the next time step. */
  void solve(double tau, double t);

  /*! @brief Assembles the left-hand-side matrix for the given time step. */
  void assemble_matrix(double tau);

  /*! @brief Assembles the right-hand-side vectors for the given time step at the given time. */
  void assemble_rhs(double tau, double t);

  /*! @brief Assembles matrix and right-hand-side. */
  void assemble(double tau, double t);

  static double get_C(const Edge &e);

  static double get_L(const Edge &e);

  static double get_R(const Edge &e);

  /*! @brief Returns the current p and q values at a given vertex. */
  void get_1d_pq_values_at_vertex(const Vertex &v, double &p, double &q) const;

  /*! @brief Evaluates p and q of the current solution on the edge e parametrized on [0, 1] at \f$ s \in [0,1] \f$. */
  void evaluate_1d_pq_values(const Edge &e, double s, double &p, double &q) const;

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
