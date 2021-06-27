////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include "mpi.h"
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"

#include "create_3_vessel_network.hpp"

namespace mc = macrocirculation;

/*! @brief Checks if we get the same solution values for A and Q as always.
 *         This should prevent unintended changes to the solver due to refactorings.
 *         Of course the resulting values should not be taken too seriously and might change if bugs are found.
 */
TEST_CASE("NonlinearSolverBitwise", "[NonlinearSolverBitwise]") {
  const double t_end = 1.2;
  const std::size_t max_iter = 160000000;
  const size_t degree = 2;
  const double tau = 5e-5;

  // create_for_node the ascending aorta
  auto graph = test_macrocirculation::util::create_3_vessel_network();

  // partition graph
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // initialize dof map
  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> solver(MPI_COMM_WORLD, graph, dof_map);
  solver.set_tau(tau);
  solver.use_ssp_method();

  for (std::size_t it = 0; it < max_iter; it += 1) {
    solver.solve();

    if (solver.get_time() > t_end + 1e-12)
      break;
  }

  // values of A we have calculated previously:
  const std::vector<double> edge_id_to_A{
    6.5155018925380483e+00,
    6.1024958654010675e+00,
    1.7777985228282513e+00};

  // values of Q we have calculated previously:
  const std::vector<double> edge_id_to_Q{
    4.1595860742808094e+02,
    3.1718725628944367e+02,
    9.1211576554297366e+01};

  for (auto e_id : graph->get_active_edge_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    auto &edge = *graph->get_edge(e_id);
    double A, Q;
    solver.evaluate_1d_AQ_values(edge, 0.5, A, Q);
    REQUIRE(A == Approx(A).epsilon(1e-16));
    REQUIRE(Q == Approx(Q).epsilon(1e-16));
  }
}