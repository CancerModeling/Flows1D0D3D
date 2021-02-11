////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
#define TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP

#include <vector>
#include <memory>
#include <functional>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class VesselDataStorage;
class DofMapNetwork;

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
  RightHandSideEvaluator(std::shared_ptr<GraphStorage> graph, std::shared_ptr<VesselDataStorage> vessel_data, std::shared_ptr<DofMapNetwork> dof_map);

  void evaluate(double t, const std::vector<double> &u_prev, std::vector<double> &rhs);

private:
  /*! @brief The current domain for solving the equation. */
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief Stores the physical data for a given vessel. */
  std::shared_ptr<VesselDataStorage> d_vessel_data;

  /*! @brief The dof map for our domain */
  std::shared_ptr<DofMapNetwork> d_dof_map;

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

  // All the parameters, which are now kept constant, but should change for vessel types.
  // TODO: These values should be provided by the network.
  double d_mu;
  double d_gamma;

  /*! @brief Recalculates the current fluxes from the previous time step.
   *
   * @param t       The current time for the inflow boundary conditions.
   * @param u_prev  The solution for which we calculate the fluxes.
   */
  void calculate_fluxes(double t, const std::vector<double> &u_prev);

  /*! @brief Assembles from the fluxes and the previous values a new right hand side function. */
  void calculate_rhs(const std::vector<double> &u_prev, std::vector<double> &rhs);

  /*! @brief Applies the inverse mass to our right hand side and saves the result to u_now. */
  void apply_inverse_mass(std::vector<double> &rhs);
};

}

#endif //TUMORMODELS_RIGHT_HAND_SIDE_EVALUATOR_HPP
