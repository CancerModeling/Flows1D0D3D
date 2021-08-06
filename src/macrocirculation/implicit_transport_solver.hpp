////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_IMPLICIT_TRANSPORT_SOLVER_HPP
#define TUMORMODELS_IMPLICIT_TRANSPORT_SOLVER_HPP

#include "mpi.h"

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class DofMap;
class PetscVec;
class PetscMat;
class PetscKsp;
class Vertex;
class Edge;
class QuadratureFormula;

class UpwindProvider {
public:
  virtual ~UpwindProvider() = default;

  virtual void init()

  /*! @brief Returns Q and A evaluated at the quadrature points of a micro-edge. */
  virtual void get_values_at_qp(double t,
                                const Edge &edge,
                                size_t micro_edge,
                                const QuadratureFormula &fe,
                                std::vector<double> &v_qp) = 0;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  virtual void get_upwinded_values(double t,
                                   const Edge &edge,
                                   std::vector<double> &v_qp) = 0;

  virtual void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) = 0;
};


class ImplicitTransportSolver {
public:
  ImplicitTransportSolver(MPI_Comm comm,
                          std::shared_ptr<GraphStorage> graph,
                          std::shared_ptr<DofMap> dof_map,
                          std::shared_ptr<UpwindProvider> upwind_provider,
                          size_t degree);

  /*! @brief Assembles matrix and right-hand-side. */
  void assemble(double tau, double t);

  void solve(double tau, double t);

  const PetscVec &get_solution() const { return *u; }

private:
  /*! @brief Assembles the left-hand-side matrix for the given time step. */
  void assemble_matrix(double tau, double t);

  /*! @brief Assembles the right-hand-side vectors for the given time step at the given time. */
  void assemble_rhs(double tau, double t);

  void assemble_matrix_cells(double tau, double t);

  void assemble_rhs_cells(double tau, double t);

  void assemble_rhs_inflow(double tau, double t);

  void assemble_matrix_outflow(double tau, double t);

  void assemble_matrix_inner_boundaries(double tau, double t);

private:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<DofMap> d_dof_map;

  std::shared_ptr<UpwindProvider> d_upwind_provider;

  size_t d_degree;

  std::shared_ptr<PetscVec> u;
  std::shared_ptr<PetscVec> rhs;
  std::shared_ptr<PetscMat> A;

  std::shared_ptr<PetscVec> mass;

  std::shared_ptr<PetscKsp> linear_solver;

  std::function<double(double)> inflow_function = [](double t) { return 0.5*std::sin(M_PI * t); };
};

} // namespace macrocirculation

#endif //TUMORMODELS_IMPLICIT_TRANSPORT_SOLVER_HPP
