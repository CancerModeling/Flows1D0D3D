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

class NonlinearLinearCoupling {
public:
  NonlinearLinearCoupling(
    MPI_Comm comm,
    std::shared_ptr<GraphStorage> graph_nl,
    std::shared_ptr<GraphStorage> graph_li);

  void add_coupled_vertices(const std::string &name);

  void add_coupled_vertices(const std::string &name_nl, const std::string &name_li);

  void update_linear_solver( const ExplicitNonlinearFlowSolver& nonlinear_solver, ImplicitLinearFlowSolver& linear_solver );

  void update_nonlinear_solver( const ImplicitLinearFlowSolver& linear_solver,  ExplicitNonlinearFlowSolver& nonlinear_solver );

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
