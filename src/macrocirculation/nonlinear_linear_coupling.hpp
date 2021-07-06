////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_NONLINEAR_LINEAR_COUPLING_HPP
#define TUMORMODELS_NONLINEAR_LINEAR_COUPLING_HPP

#include <memory>

#include "communication/buffer.hpp"

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class ExplicitNonlinearFlowSolver;
class ImplicitLinearFlowSolver;
class Vertex;

/*! @brief Couples a graph using a nonlinear explicit flow solver with a graph based on a linear implicit flow solver by updating the respective characteristic boundary conditions.
 *
 * For this coupling we use the linear boundary values (p,q) for a nonlinear characteristic boundary condition on (A,Q)
 * and the nonlinear boundary values (A,Q) for a linear characteristic boundary condition on (p,q).
 */
class NonlinearLinearCoupling {
public:
  /*!
   * @param comm The MPI communicator.
   * @param graph_nl The graph for the nonlinear explicit flow model.
   * @param graph_li The graph for the linear implicit flow model.
   */
  NonlinearLinearCoupling(
    MPI_Comm comm,
    std::shared_ptr<GraphStorage> graph_nl,
    std::shared_ptr<GraphStorage> graph_li);

  /*! @brief Couples vertices from the graph with nonlinear flow to the graph with linear flow.
   *
   * @param name A vertex name appearing in both graphs, which should be coupled.
   */
  void add_coupled_vertices(const std::string &name);

  /*! @brief Couples vertices from the graph with nonlinear flow to the graph with linear flow.
   *
   * @param name_nl The vertex name on the nonlinear flow graph.
   * @param name_li The vertex name on the linear flow graph.
   */
  void add_coupled_vertices(const std::string &name_nl, const std::string &name_li);

  /*! @brief Updates the (characteristic) boundary conditions of the linear solver with boundary data from the nonlinear solver.
   *
   * @param nonlinear_solver The nonlinear explicit solver.
   * @param linear_solver The linear implicit solver.
   */
  void update_linear_solver(const ExplicitNonlinearFlowSolver &nonlinear_solver, ImplicitLinearFlowSolver &linear_solver);

  /*! @brief Updates the (characteristic) boundary conditions of the nonlinear solver with boundary data from the linear solver.
   *
   * @param linear_solver The linear implicit solver.
   * @param nonlinear_solver The nonlinear explicit solver.
   */
  void update_nonlinear_solver(const ImplicitLinearFlowSolver &linear_solver, ExplicitNonlinearFlowSolver &nonlinear_solver);

private:
  static int get_rank(GraphStorage &graph, Vertex &v);

  struct CoupledVertices {
    size_t vertex_id_1;
    size_t vertex_id_2;
  };

  std::vector<CoupledVertices> coupled_vertices;

private:
  std::shared_ptr<GraphStorage> d_graph_nl;
  std::shared_ptr<GraphStorage> d_graph_li;

  MPI_Comm d_comm;

  BufferSystem d_buffer_system;
};

} // namespace macrocirculation

#endif //TUMORMODELS_NONLINEAR_LINEAR_COUPLING_HPP
