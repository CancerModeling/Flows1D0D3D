////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
#define TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP

#include "nonlinear_flow_upwind_evaluator.hpp"
#include <functional>
#include <memory>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class Vertex;
class Edge;

/*! @brief Functional to evaluate the right-hand-side S. */
class default_S {
public:
  explicit default_S(double phi);

  /*! @brief Evaluates the Q- and A-component of the right-hand side S,
   *         given Q and A at eacht of the quadrature points.
   */
  void operator()(double,
                  const Edge &,
                  const std::vector<double> &,
                  const std::vector<double> &Q,
                  const std::vector<double> &A,
                  std::vector<double> &S_Q_out,
                  std::vector<double> &S_A_out) const;

private:
  /*! @brief Wall permeability [cm^2 s^{-1}]. */
  double d_phi;
};

/*! @brief Assembles the inverse mass. WARNING: Assumes legendre basis! */
void assemble_inverse_mass(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, std::vector<double> &inv_mass);

/*! @brief Our flow equation d/dt u + d/dz F(u) = 0, can be written with a DG ansatz as,
 *              d/dt u - M^{-1}( (F(u), d/dz v) - F(u(x_r))v(x_r) + F(u(x_l))v(x_l) )  = 0,
 *         where M^{-1} is the inverse of the mass matrix, v a test function and x_r and x_l
 *         the boundaries of an interval.
 *         This class evaluates
 *              M^{-1}( (F(u), d/dz v) - F(u(x_r))v(x_r) + F(u(x_l))v(x_l) )
 *         and can be easily used in a time integrator.
 */
class RightHandSideEvaluator {
public:
  RightHandSideEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map, std::size_t degree);

  void evaluate(double t, const std::vector<double> &u_prev, std::vector<double> &rhs);

  /*! @brief Function type to evaluate a vectorial quantity at all the quadrature points in one go:
   *         - the 1st argument is the current time,
   *         - the 2nd argument the vessel,
   *         - the 3rd argument are the quadrature points on the edge,
   *         - the 4th and 5th components are for the flow Q and the area A, while
   *         - the 5th and 6th are the components of the vector evaluated at the quadrature points.
   */
  using VectorEvaluator = std::function<void(double,
                                             const Edge &,
                                             const std::vector<double> &,
                                             const std::vector<double> &,
                                             const std::vector<double> &,
                                             std::vector<double> &,
                                             std::vector<double> &)>;

  /*! @brief Sets the right-hand side S. */
  void set_rhs_S(VectorEvaluator S_evaluator);

private:
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMap> d_dof_map;

  NonlinearFlowUpwindEvaluator d_flow_upwind_evaluator;

  /*! @brief Evaluates the Q- and A-component of the right-hand side S,
   *         given Q and A at eacht of the quadrature points.
   */
  VectorEvaluator d_S_evaluator;

  /*! @brief The degree of the finite-element shape functions. */
  std::size_t d_degree;

  /*! @brief Our current inverse mass vector, defining the diagonal inverse mass matrix. */
  std::vector<double> d_inverse_mass;

  /*! @brief Assembles from the fluxes and the previous values a new right hand side function. */
  void calculate_rhs(double t, const std::vector<double> &u_prev, std::vector<double> &rhs);

  /*! @brief Applies the inverse mass to our right hand side and saves the result to u_now. */
  void apply_inverse_mass(std::vector<double> &rhs);
};

} // namespace macrocirculation

#endif //TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
