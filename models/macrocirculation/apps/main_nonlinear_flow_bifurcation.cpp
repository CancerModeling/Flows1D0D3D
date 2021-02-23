////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <memory>
#include <chrono>

#include "macrocirculation/graph_data_writer.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/vessel_data_storage.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/communicator.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/graph_partitioner.hpp"


namespace lm = libMesh;
namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  lm::LibMeshInit init(argc, argv);

  const double t_end = 4e-1;
  const std::size_t max_iter = 160000000;

  const double tau = 2.5e-4/32/4;
  const double tau_out = 1e-3;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  const std::size_t num_edges_per_segment = 80;

  // we create data for celiac ii
  auto vessel_data = std::make_shared< mc::VesselDataStorage > ();
  std::size_t main_vessel_id = vessel_data->add_parameter({
    1706.7e2, // 1706.7 hPa
    0.13, // 0.13 cm^2
    1.028,   // 1.028 kg/cm^3,
  });

  std::size_t other_vessel_id = vessel_data->add_parameter({
    1706.7e2/std::sqrt(2),
    0.13/2,
  1.028,
  });

  const double scale = 1;

  // create the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  // create vertices:
  auto start = graph->create_vertex(lm::Point(0, 0, 0));
  auto midpoint1 = graph->create_vertex(lm::Point(1*scale, 0, 0));
  auto upper_point = graph->create_vertex(lm::Point(1.5*scale, +0.5*scale, 0));
  auto lower_point= graph->create_vertex(lm::Point(1.5*scale, -0.5*scale, 0));
  auto midpoint2 = graph->create_vertex(lm::Point(2*scale, 0, 0));
  auto endpoint = graph->create_vertex(lm::Point(3*scale, 0, 0));
  // connect vertices:
  graph->line_to(*start, *midpoint1, main_vessel_id, num_edges_per_segment);
  graph->line_to(*midpoint1, *upper_point, other_vessel_id, num_edges_per_segment);
  graph->line_to(*midpoint1, *lower_point, other_vessel_id, num_edges_per_segment);
  graph->line_to(*lower_point, *midpoint2, other_vessel_id, num_edges_per_segment);
  graph->line_to(*upper_point, *midpoint2, other_vessel_id, num_edges_per_segment);
  graph->line_to(*midpoint2, *endpoint, main_vessel_id, num_edges_per_segment);
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // set inflow boundary conditions
  start->set_to_inflow(mc::heart_beat_inflow());

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> solver(MPI_COMM_WORLD, graph, vessel_data);
  solver.set_tau(tau);
  solver.use_ssp_method();

  std::vector<double> Q_vertex_values(graph->num_edges() * 2, 0);
  std::vector<double> A_vertex_values(graph->num_edges() * 2, 0);
  std::vector<double> p_vertex_values(graph->num_edges() * 2, 0);
  std::vector<double> p_static_vertex_values(graph->num_edges() * 2, 0);

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {

    solver.solve();

    if (it % output_interval == 0) {
      std::cout << "iter " << it << std::endl;

      solver.get_communicator().gather(0, solver.get_solution());

      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
      {
        // save solution
        solver.get_solution_on_vertices(Q_vertex_values, A_vertex_values);
        solver.get_total_pressure_on_vertices(p_vertex_values);
        solver.get_static_pressure_on_vertices(p_static_vertex_values);
        mc::GraphDataWriter writer;
        writer.add_vertex_data("Q", Q_vertex_values);
        writer.add_vertex_data("A", A_vertex_values);
        writer.add_vertex_data("p_total", p_vertex_values);
        writer.add_vertex_data("p_static", p_static_vertex_values);
        writer.write_vtk("bifurcation_solution", *graph, it);
      }
    }

    // break
    if (solver.get_time() > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms =  std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout  << "time = " << elapsed_ms*1e-6 << " s" << std::endl;
}
