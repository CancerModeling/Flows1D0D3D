////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
#define TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP

#include <vector>
#include <memory>
#include <functional>
#include <mpi.h>

#include "libmesh/point.h"

namespace macrocirculation {

using libMesh::Point;

// forward declarations
class GraphStorage;
class VesselDataStorage;
class DofMapNetwork;
class Communicator;

/*! @brief Functional to evaluate the right-hand-side S. */
class default_S {
public:
  default_S(double mu, double gamma, double phi);

  /*! @brief Evaluates the Q- and A-component of the right-hand side S,
   *         given Q and A at eacht of the quadrature points.
   */
  void operator()(double, const std::vector<Point>&, const std::vector<double>& Q, const std::vector<double>& A, std::vector<double>& S_Q_out, std::vector<double>& S_A_out) const;

private:
  /*! @brief Blood viscosity [m Pa s]. */
  double d_mu;

  /*! @brief Shape of velocity profile. */
  double d_gamma;

  /*! @brief Wall permeability [cm^2 s^{-1}]. */
  double d_phi;
};

/*! @brief Assembles the inverse mass. WARNING: Assumes legendre basis! */
template <std::size_t degree>
void assemble_inverse_mass(const GraphStorage &graph, const DofMapNetwork &dof_map, std::vector<double> &inv_mass);

/*! @brief Our flow equation d/dt u + d/dz F(u) = 0, can be written with a DG ansatz as,
 *              d/dt u - M^{-1}( (F(u), d/dz v) - F(u(x_r))v(x_r) + F(u(x_l))v(x_l) )  = 0,
 *         where M^{-1} is the inverse of the mass matrix, v a test function and x_r and x_l
 *         the boundaries of an interval.
 *         This class evaluates
 *              M^{-1}( (F(u), d/dz v) - F(u(x_r))v(x_r) + F(u(x_l))v(x_l) )
 *         and can be easily used in a time integrator.
 */
template <std::size_t degree>
class RightHandSideEvaluator {
public:
  RightHandSideEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data, std::shared_ptr<DofMapNetwork> dof_map);

  void evaluate(double t, std::vector<double> &u_prev, std::vector<double> &rhs);

  /*! @brief Function type to evaluate a vectorial quantity at all the quadrature points in one go:
   *         - the 1st argument is the current time,
   *         - the 2nd argument are the quadrature points on the domain,
   *         - the 3nd and 4rd components are for the flow Q and the area A, while
   *         - the 5th and 6th are the components of the vector evaluated at the quadrature points.
   */
  using VectorEvaluator = std::function< void(double, const std::vector<Point>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&) >;

  /*! @brief Sets the right-hand side S. */
  void set_rhs_S(VectorEvaluator S_evaluator);

private:
  /*! @brief The mpi communicator. */
  MPI_Comm d_comm;

  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief Stores the physical data for a given vessel. */
  std::shared_ptr<VesselDataStorage> d_vessel_data;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMapNetwork> d_dof_map;

  /*! @brief The mpi-communication routine to update the ghost layers. */
  std::unique_ptr<Communicator> d_communicator;

  /*! @brief Evaluates the Q- and A-component of the right-hand side S,
   *         given Q and A at eacht of the quadrature points.
   */
  VectorEvaluator d_S_evaluator;

  /*! @brief The upwinded Q-flow from the left boundary of a macroedge,
   *         ordered such that d_Q_up_el[edge_id] contains the respective flux.
   */
  std::vector<double> d_Q_up_el;

  /*! @brief The upwinded Q-flow from the right boundary of a macroedge,
   *         ordered such that d_Q_up_er[edge_id] contains the respective flux.
   */
  std::vector<double> d_Q_up_er;

  /*! @brief The upwinded A-flow from the left boundary of a macroedge,
   *         ordered such that d_A_up_el[edge_id] contains the respective flux.
   */
  std::vector<double> d_A_up_el;

  /*! @brief The upwinded A-flow from the right boundary of a macroedge,
   *         ordered such that d_Q_up_er[edge_id] contains the respective flux.
   */
  std::vector<double> d_A_up_er;

  /*! @brief Our current inverse mass vector, defining the diagonal inverse mass matrix. */
  std::vector<double> d_inverse_mass;

  /*! @brief We cache the active vertex ids. */
  std::vector<std::size_t> d_active_vertex_ids;

  /*! @brief We cache the active edge ids. */
  std::vector<std::size_t> d_active_edge_ids;

  /*! @brief Recalculates the current fluxes from the previous time step.
   *
   * @param t       The current time for the inflow boundary conditions.
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_fluxes(double t, const std::vector<double> &u_prev);

  /*! @brief Assembles from the fluxes and the previous values a new right hand side function. */
  void calculate_rhs(double t, const std::vector<double> &u_prev, std::vector<double> &rhs);

  /*! @brief Applies the inverse mass to our right hand side and saves the result to u_now. */
  void apply_inverse_mass(std::vector<double> &rhs);
};

}

#endif //TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
