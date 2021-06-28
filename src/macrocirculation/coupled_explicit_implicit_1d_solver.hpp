////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_COUPLED_EXPLICIT_IMPLICIT_1D_SOLVER_HPP
#define TUMORMODELS_COUPLED_EXPLICIT_IMPLICIT_1D_SOLVER_HPP

#include "mpi.h"
#include <memory>

namespace macrocirculation {

// forward declarations
class GraphStorage;
class DofMap;
class ExplicitNonlinearFlowSolver;
class ImplicitLinearFlowSolver;
class NonlinearLinearCoupling;

class CoupledExplicitImplicit1DSolver {
public:
  CoupledExplicitImplicit1DSolver(
    MPI_Comm comm,
    std::shared_ptr< NonlinearLinearCoupling > coupling,
    std::shared_ptr<GraphStorage> graph_nl,
    std::shared_ptr<GraphStorage> graph_li,
    size_t degree_nl,
    size_t degree_li);

  void setup(double tau);

  void solve(double tau, double t);

  std::shared_ptr<DofMap> get_explicit_dof_map() { return d_dof_map_nl; }
  std::shared_ptr<DofMap> get_implicit_dof_map() { return d_dof_map_li; }

  std::shared_ptr<ExplicitNonlinearFlowSolver> get_explicit_solver() { return d_explicit_solver; }
  std::shared_ptr<ImplicitLinearFlowSolver> get_implicit_solver() { return d_implicit_solver; }

  std::shared_ptr<NonlinearLinearCoupling> get_coupling() { return d_coupling; }

protected:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph_nl;
  std::shared_ptr<GraphStorage> d_graph_li;

  std::shared_ptr<DofMap> d_dof_map_nl;
  std::shared_ptr<DofMap> d_dof_map_li;

  std::shared_ptr<ExplicitNonlinearFlowSolver> d_explicit_solver;
  std::shared_ptr<ImplicitLinearFlowSolver> d_implicit_solver;

  std::shared_ptr<NonlinearLinearCoupling> d_coupling;

  double d_tau;

  bool d_is_setup;
};

} // namespace macrocirculation

#endif //TUMORMODELS_COUPLED_EXPLICIT_IMPLICIT_1D_SOLVER_HPP
