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

  /*! @brief Currently our geometry consists of 3 vessels, s.t.
   *         +--------+-[Aneurysm]-+-------+
   *                  in          out
   *         This function allows us to identify, if we want the data at in or at out.
   */
  enum class Outlet { in,
                      out };

  enum class Inlet { in,
                     out };

  /*! @brief Propagates the solver and its solution one time step further. */
  void solve();

  /*! @brief Returns the coupling data at the given outlet (either in or out, see above). */
  Result get_result(Outlet outlet);

  /*! @brief Returns the outer values just to check. */
  Result get_result_outer(Outlet outlet);

  void set_result(Outlet outlet, double p, double q);

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

  std::shared_ptr<Vertex> v_coupling_1_outer;
  std::shared_ptr<Vertex> v_coupling_1_inner;
  std::shared_ptr<Vertex> v_coupling_2_inner;
  std::shared_ptr<Vertex> v_coupling_2_outer;

  std::shared_ptr<Edge> edge0;
  std::shared_ptr<Edge> edge1;

  std::shared_ptr<GraphPVDWriter> pvd_writer;
  std::shared_ptr<GraphCSVWriter> csv_writer;
};

} // namespace macrocirculation

#endif