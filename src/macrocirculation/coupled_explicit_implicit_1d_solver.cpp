//
// Created by wagneran on 28.06.21.
//

#include "coupled_explicit_implicit_1d_solver.hpp"

#include "petsc.h"
#include <cmath>
#include <utility>

#include "dof_map.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "nonlinear_linear_coupling.hpp"

namespace macrocirculation {

CoupledExplicitImplicit1DSolver::CoupledExplicitImplicit1DSolver(
  MPI_Comm comm,
  std::shared_ptr<NonlinearLinearCoupling> coupling,
  std::shared_ptr<GraphStorage> graph_nl,
  std::shared_ptr<GraphStorage> graph_li,
  size_t degree_nl,
  size_t degree_li)
    : d_comm(comm),
      d_graph_nl(std::move(graph_nl)),
      d_graph_li(std::move(graph_li)),
      d_dof_map_nl(nullptr),
      d_dof_map_li(nullptr),
      d_explicit_solver(nullptr),
      d_implicit_solver(nullptr),
      d_coupling(std::move(coupling)),
      d_tau(NAN),
      d_is_setup(false) {

  d_graph_nl->finalize_bcs();
  d_graph_li->finalize_bcs();

  d_dof_map_nl = std::make_shared<DofMap>(d_graph_nl->num_vertices(), d_graph_nl->num_edges());
  d_dof_map_nl->create(comm, *d_graph_nl, 2, degree_nl, false);

  d_dof_map_li = std::make_shared<DofMap>(d_graph_li->num_vertices(), d_graph_li->num_edges());
  d_dof_map_li->create(comm, *d_graph_li, 2, degree_li, true);

  d_explicit_solver = std::make_shared<ExplicitNonlinearFlowSolver>(MPI_COMM_WORLD, d_graph_nl, d_dof_map_nl, degree_nl);
  d_implicit_solver = std::make_shared<ImplicitLinearFlowSolver>(MPI_COMM_WORLD, d_graph_li, d_dof_map_li, degree_li);

  d_explicit_solver->use_explicit_euler_method();
}

void CoupledExplicitImplicit1DSolver::setup(double tau) {
  d_explicit_solver->set_tau(tau);
  d_implicit_solver->setup(tau);
  d_tau = tau;
  d_is_setup = true;
}

void CoupledExplicitImplicit1DSolver::solve(double tau, double t) {
  if (!d_is_setup || std::abs(tau - d_tau) > 1e-16)
    throw std::runtime_error("missing call to setup routines");

  d_explicit_solver->solve();
  d_coupling->update_linear_solver(*d_explicit_solver, *d_implicit_solver);
  d_implicit_solver->solve(d_tau, t + d_tau);
  d_coupling->update_nonlinear_solver(*d_implicit_solver, *d_explicit_solver);
}

} // namespace macrocirculation
