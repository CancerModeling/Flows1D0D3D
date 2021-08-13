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
#include <random>

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
class FlowUpwindEvaluator;
class ExplicitNonlinearFlowSolver;

class UpwindProvider {
public:
  virtual ~UpwindProvider() = default;

  virtual void init(double t, const std::vector<double> &u) {}
  virtual void init(double t, const PetscVec &u) {}

  /*! @brief Returns Q and A evaluated at the quadrature points of a micro-edge. */
  virtual void get_values_at_qp(double t,
                                const Edge &edge,
                                size_t micro_edge,
                                const QuadratureFormula &fe,
                                std::vector<double> &v_qp) const = 0;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  virtual void get_upwinded_values(double t,
                                   const Edge &edge,
                                   std::vector<double> &v_qp) const = 0;

  virtual void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const = 0;
};

class ConstantUpwindProvider : public UpwindProvider {
public:
  explicit ConstantUpwindProvider(double speed);

  ~ConstantUpwindProvider() override;

  void get_values_at_qp(double t,
                        const Edge &edge,
                        size_t micro_edge,
                        const QuadratureFormula &qf,
                        std::vector<double> &v_qp) const override;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const override;

  void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const override;

private:
  double d_speed;
};

class UpwindProviderNonlinearFlow : public UpwindProvider {
public:
  explicit UpwindProviderNonlinearFlow(std::shared_ptr<FlowUpwindEvaluator> evaluator, std::shared_ptr<ExplicitNonlinearFlowSolver> solver);

  ~UpwindProviderNonlinearFlow() override = default;

  void init(double t, const std::vector<double> &u) override;

  void get_values_at_qp(double t,
                        const Edge &edge,
                        size_t micro_edge,
                        const QuadratureFormula &qf,
                        std::vector<double> &v_qp) const override;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const override;

  void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const override;

private:
  std::shared_ptr<FlowUpwindEvaluator> d_evaluator;
  std::shared_ptr<ExplicitNonlinearFlowSolver> d_solver;
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

  /*! @brief Solves the system. For the _next_ time step t.
   *
   * @param tau The current time step width.
   * @param t The (future) time for the next time step. */
  void solve(double tau, double t);

  const PetscVec &get_solution() const { return *u; }

private:
  /*! @brief Assembles the left-hand-side matrix for the given time step. */
  void assemble_matrix(double tau, double t, const UpwindProvider& upwind_provider);

  /*! @brief Assembles the right-hand-side vectors for the given time step at the given time. */
  void assemble_rhs(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_matrix_cells(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_rhs_cells(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_rhs_inflow(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_matrix_outflow(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_matrix_inner_boundaries(double tau, double t, const UpwindProvider& upwind_provider);

  void assemble_matrix_nfurcations(double tau, double t, const UpwindProvider& upwind_provider);

  /*! @brief Assembles the matrix with different upwindings, so that the sparsity pattern does not change,
   *         when the velocity field changes.  */
  void sparsity_pattern();

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

  // std::function<double(double)> inflow_function = [](double t) { return 0.5*std::sin(M_PI * t); };
  std::function<double(double)> inflow_function = [](double t) {
    const double delta = 0.05;
    if (t < delta)
      return -2 * std::pow(t / delta, 3) + 3 * std::pow(t / delta, 2);
    else
      return 1.;
  };
};

} // namespace macrocirculation

#endif //TUMORMODELS_IMPLICIT_TRANSPORT_SOLVER_HPP
