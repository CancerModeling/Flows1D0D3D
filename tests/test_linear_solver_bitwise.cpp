////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include "mpi.h"
#include <iomanip>
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"

#include "create_3_vessel_network.hpp"

namespace mc = macrocirculation;

/*! @brief Checks if we get the same solution values for p and q as always.
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
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, true);

  // configure solver
  mc::LinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map, degree);
  solver.get_solver().set_pc_to_jacobi();
  solver.setup(tau);

  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {
    t += tau;
    solver.solve(tau, it * tau);

    if (t > t_end + 1e-12)
      break;
  }

  // values of A we have calculated previously:
  const std::vector<double> edge_id_to_p{
    1.4862924418636430e+02,
    1.4995440475551933e+02,
    1.5002923944324172e+02};

  // values of Q we have calculated previously:
  const std::vector<double> edge_id_to_q{
    4.1768385246597893e+02,
    3.1987614522229723e+02,
    9.2213401604731729e+01};

  for (auto e_id : graph->get_active_edge_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    auto &edge = *graph->get_edge(e_id);
    double p, q;
    solver.evaluate_1d_pq_values(edge, 0.5, p, q);
    std::cout << e_id << std::scientific << std::setprecision(16) << " p = " << p << " q = " << q << std::endl;
    REQUIRE(p == Approx(edge_id_to_p[e_id]).epsilon(1e-10));
    REQUIRE(q == Approx(edge_id_to_q[e_id]).epsilon(1e-10));
  }
}