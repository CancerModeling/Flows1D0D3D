////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_IMPLICIT_ADVECTION_SOLVER_H
#define TUMORMODELS_IMPLICIT_ADVECTION_SOLVER_H

#include <functional>
#include <memory>
#include <mpi.h>
#include <utility>

namespace macrocirculation {

class GraphStorage;

class ImplicitAdvectionSolver {
public:
  /*! @brief Type of the function providing inflow boundary values.
   *         The argument is the current time.
   */
  using InflowValueFct = std::function<double(double)>;

  ImplicitAdvectionSolver(std::shared_ptr<GraphStorage> graph, std::size_t degree, double tau, double t_end, double velocity, InflowValueFct inflow_value_fct);

  void solve() const;

  void set_output_interval(std::size_t interval) { d_output_interval = interval; }

private:
  MPI_Comm d_comm;

  std::size_t d_degree;

  double d_tau;
  double d_t_end;
  double d_velocity;
  std::size_t d_output_interval;
  InflowValueFct d_inflow_value_fct;

  std::shared_ptr<GraphStorage> d_graph;
};

} // namespace macrocirculation

#endif //TUMORMODELS_IMPLICIT_ADVECTION_SOLVER_H
