////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_ADVECTION_SOLVER_HPP
#define TUMORMODELS_ADVECTION_SOLVER_HPP

#include <functional>
#include <utility>
#include <memory>

namespace macrocirculation {

class GraphStorage;

template < std::size_t DEGREE >
class AdvectionSolver {
public:
  /*! @brief Type of the function providing inflow boundary values.
   *         The argument is the current time.
   */
  using InflowValueFct = std::function<double(double)>;

  AdvectionSolver(std::shared_ptr<GraphStorage> graph, double tau, double t_end, double velocity, InflowValueFct inflow_value_fct);

  void solve() const;

  void set_output_interval(std::size_t interval){ d_output_interval = interval; }

private:
  double d_tau;
  double d_t_end;
  double d_velocity;
  std::size_t d_output_interval;
  InflowValueFct d_inflow_value_fct;

  std::shared_ptr<GraphStorage> d_graph;
};

} // namespace macrocirculation

#endif //TUMORMODELS_ADVECTION_SOLVER_HPP
