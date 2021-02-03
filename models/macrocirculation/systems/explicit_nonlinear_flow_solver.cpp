////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "explicit_nonlinear_flow_solver.h"

namespace macrocirculation {

constexpr std::size_t degree = 2;

ExplicitNonlinearFlowSolver::ExplicitNonlinearFlowSolver(std::shared_ptr<GraphStorage> graph)
    : d_graph(graph) {
  // TODO: Initialize the rest
}

void ExplicitNonlinearFlowSolver::solve() {
  // TODO: Implement
}

void ExplicitNonlinearFlowSolver::calculate_fluxes() {
  // TODO: Implement
}

} // namespace macrocirculation