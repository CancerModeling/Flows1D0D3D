////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_LINEARIZED_SOLVER_H
#define SIMPLE_LINEARIZED_SOLVER_H

#include <memory>
#include <vector>
#include <functional>
#include <mpi.h>

namespace macrocirculation {

class GraphStorage;
class ImplicitLinearFlowSolver;
class Edge;
class Vertex;
class DofMap;
class GraphPVDWriter;
class GraphCSVWriter;

/*! This is a simplified interface for calling our code from an LBM aneurysm solver. */
class SimpleLinearizedSolver {
public:
  SimpleLinearizedSolver(MPI_Comm comm, const std::string & filepath, const std::string& folder, const std::string& name, double tau);

  SimpleLinearizedSolver(const std::string & filepath, const std::string& folder, const std::string& name, double tau = 1e-5);

  struct Result {
    /*! @brief Vessel area in [cm^2]. */
    double a;
    /*! @brief Pressure in [dyne / cm^2]. */
    double p;
    /*! Flow in [cm^3 / s]. */
    double q;
  };

  /*! @brief Propagates the solver and its solution one time step further. */
  void solve();

  /*! @brief Returns the coupling data at the given outlet specified with an index. */
  Result get_result(int outlet);

  std::vector<std::array<double, 3>> get_points();

  /*! @brief Returns the outer values just to check. */
  Result get_result_outer(int outlet);

  void set_result(int outlet, double p, double q);

  /*! @brief Writes the solution to the disk for debugging purposes.
   *         Output is pretty slows, so don't call this too often.
   */
  void write();

  /*! @brief Sets the time step size. */
  void set_tau(double tau);

  /*! @brief Sets the inflow function. */
  void set_inflow(const std::function< double(double)> & fun);

  /*! @brief Sets the outflow to an RCR model. */
  void set_outflow_rcr(const double R, const double C);

  /*! @brief Returns the number of coupling points found in the 1D geometry. */
  size_t get_num_coupling_points() const;

  /*! @brief Returns the simulation time, to which the currently saved solution and the get_result outputs correspond. */
  double get_current_t() const;

  void set_t(double t);

private:
  Result get_result(const Vertex &vertex, const Edge &edge);

  void setup();

private:
  MPI_Comm comm;

  double tau{};
  double t{};

  std::shared_ptr<GraphStorage> graph;
  std::shared_ptr<DofMap> dof_map;
  std::shared_ptr<ImplicitLinearFlowSolver> solver;

  std::vector<std::shared_ptr<Vertex>> v_coupling_outer;
  std::vector<std::shared_ptr<Vertex>> v_coupling_inner;

  std::vector<std::shared_ptr<Edge>> edge;

  std::shared_ptr<GraphPVDWriter> pvd_writer;
  std::shared_ptr<GraphCSVWriter> csv_writer;

  size_t num_coupling_points;

  double rescale_q(const Vertex &vertex, const Edge &edge, double q);
};

} // namespace macrocirculation

#endif