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
#include <random>
#include <vector>

#include "petsc_assembly_blocks.hpp"

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
class NonlinearFlowUpwindEvaluator;
class LinearizedFlowUpwindEvaluator;
class ExplicitNonlinearFlowSolver;
class ImplicitLinearFlowSolver;
class Point;

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

  virtual void get_0d_pressures(double t, const Vertex &v, std::vector<double> &p_c) const = 0;
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

  void get_0d_pressures(double t, const Vertex &v, std::vector<double> &p_c) const override { throw std::runtime_error("not implemented yet"); }

private:
  double d_speed;
};

class EmbeddedUpwindProvider : public UpwindProvider {
public:
  explicit EmbeddedUpwindProvider(std::shared_ptr<GraphStorage> graph, std::function<double(double, const Point &)> field, std::function<void(double, const Vertex &, std::vector<double> &p_c)> field_0d);

  ~EmbeddedUpwindProvider() override;

  void get_values_at_qp(double t,
                        const Edge &edge,
                        size_t micro_edge,
                        const QuadratureFormula &qf,
                        std::vector<double> &v_qp) const override;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const override;

  void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const override;

  void get_0d_pressures(double t, const Vertex &v, std::vector<double> &p_c) const override { d_0d_pressure_field(t, v, p_c); }

private:
  std::shared_ptr<GraphStorage> d_graph;

  std::function<double(double, const Point &)> d_field;

  std::function<void(double, const Vertex &, std::vector<double> &p_c)> d_0d_pressure_field;
};

class UpwindProviderNonlinearFlow : public UpwindProvider {
public:
  explicit UpwindProviderNonlinearFlow(std::shared_ptr<NonlinearFlowUpwindEvaluator> evaluator, std::shared_ptr<ExplicitNonlinearFlowSolver> solver);

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

  void get_0d_pressures(double t, const Vertex &v, std::vector<double> &p_c) const override;

private:
  std::shared_ptr<NonlinearFlowUpwindEvaluator> d_evaluator;
  std::shared_ptr<ExplicitNonlinearFlowSolver> d_solver;
};

class UpwindProviderLinearizedFlow : public UpwindProvider {
public:
  UpwindProviderLinearizedFlow(std::shared_ptr<GraphStorage> graph, std::shared_ptr<LinearizedFlowUpwindEvaluator> evaluator, std::shared_ptr<ImplicitLinearFlowSolver> solver);

  ~UpwindProviderLinearizedFlow() override = default;

  void init(double t, const std::vector<double> &u) override;

  void init(double t, const PetscVec &u) override;

  void get_values_at_qp(double t,
                        const Edge &edge,
                        size_t micro_edge,
                        const QuadratureFormula &qf,
                        std::vector<double> &v_qp) const override;

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const override;

  void get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const override;

  void get_0d_pressures(double t, const Vertex &v, std::vector<double> &p_c) const override;

private:
  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<LinearizedFlowUpwindEvaluator> d_evaluator;

  std::shared_ptr<ImplicitLinearFlowSolver> d_solver;
};

class ImplicitTransportSolver {
public:
  ImplicitTransportSolver(MPI_Comm comm,
                          std::shared_ptr<GraphStorage> graph,
                          std::shared_ptr<DofMap> dof_map,
                          std::shared_ptr<UpwindProvider> upwind_provider,
                          size_t degree);

  ImplicitTransportSolver(MPI_Comm comm,
                          std::vector<std::shared_ptr<GraphStorage>> graph,
                          std::vector<std::shared_ptr<DofMap>> dof_map,
                          std::vector<std::shared_ptr<UpwindProvider>> upwind_provider,
                          size_t degree);

  /*! @brief Assembles matrix and right-hand-side. */
  void assemble(double tau, double t);

  /*! @brief Solves the system. For the _next_ time step t.
   *
   * @param tau The current time step width.
   * @param t The (future) time for the next time step. */
  void solve(double tau, double t);

  const PetscVec &get_solution() const { return *u; }

  void apply_slope_limiter(double t);

  const PetscMat &get_mat() const { return *A; }

  const PetscVec &get_rhs() const { return *rhs; }

  void set_inflow_function(std::function<double(double)> inflow_function);

  // TODO: move volumes into dedicated class
  const PetscVec &get_volumes() const { return *d_volumes; }

  // TODO: move volumes into dedicated class
  const std::vector< std::shared_ptr< DofMap > > &get_dof_maps_volume() const { return d_dof_maps_volume; }

  const std::vector< std::shared_ptr< DofMap > > &get_dof_maps_transport() const { return d_dof_map; }

private:
  /*! @brief Assembles the left-hand-side matrix for the given time step. */
  void assemble_matrix(double tau, double t);

  /*! @brief Assembles the right-hand-side vectors for the given time step at the given time. */
  void assemble_rhs(double tau, double t);

  void assemble_characteristics(double tau, double t);

  void assemble_windkessel_rhs_and_matrix(double tau, double t);

  void assemble_windkessel_rhs_and_matrix(double tau,
                                          double t,
                                          const GraphStorage &graph,
                                          const DofMap &dof_map,
                                          const DofMap &dof_map_volume,
                                          const UpwindProvider &upwind_provider);

  /*! @brief Assembles the matrix with different upwindings, so that the sparsity pattern does not change,
   *         when the velocity field changes.  */
  void sparsity_pattern();

  void sparsity_pattern_characteristics();

private:
  MPI_Comm d_comm;

  std::vector<std::shared_ptr<GraphStorage>> d_graph;

  std::vector<std::shared_ptr<DofMap>> d_dof_map;

  std::vector<std::shared_ptr<UpwindProvider>> d_upwind_provider;

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

  std::vector<std::shared_ptr<DofMap>> d_dof_maps_volume;

  std::shared_ptr<PetscVec> d_volumes;

  void apply_slope_limiter(std::shared_ptr<GraphStorage> d_graph, std::shared_ptr<DofMap> d_dof_map, std::shared_ptr<UpwindProvider> upwind_provider, double t);
};

namespace implicit_transport {

void additively_assemble_rhs_cells(PetscVec &mass, PetscVec &u, PetscVec &rhs);

void additively_assemble_rhs_inflow(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, const std::function<double(double)> &inflow_function, PetscVec &rhs);

void additively_assemble_matrix(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, PetscMat &mat);

void additively_assemble_matrix_cells(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, PetscMat &mat);

void additively_assemble_matrix_outflow(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, PetscMat &mat);

void additively_assemble_matrix_inner_boundaries(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, PetscMat &mat);

void additively_assemble_matrix_nfurcations(MPI_Comm comm, double tau, double t, const UpwindProvider &upwind_provider, const DofMap &dof_map, const GraphStorage &graph, PetscMat &mat);

void assembly_kernel_nfurcation(MPI_Comm comm, double tau, const std::vector<std::vector<size_t>> &dof_indices_list, const std::vector<Edge const *> &edges, const std::vector<double> &sigma, const std::vector<double> &Q_up, const std::vector<double> &A_up, const std::vector<BoundaryPointType> &boundary_type, PetscMat &A);

} // namespace implicit_transport

} // namespace macrocirculation

#endif //TUMORMODELS_IMPLICIT_TRANSPORT_SOLVER_HPP
