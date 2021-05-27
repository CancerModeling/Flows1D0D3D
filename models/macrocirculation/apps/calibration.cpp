////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <argparse.hpp>
#include <chrono>
#include <macrocirculation/communication/mpi.hpp>
#include <memory>

#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_flow_and_concentration_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/transport.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  argparse::ArgumentParser program("Calibration script for geometries");
  program.add_argument("--input-file")
    .help("path to the input file")
    .default_value(std::string("./data/network-33-vessels.json"));
  program.add_argument("--heart-amplitude")
    .help("the amplitude of a heartbeat")
    .default_value(485.0)
    .action([](const std::string &value) { return std::stod(value); });
  program.add_argument("--tau")
    .help("time step size")
    .default_value(2.5e-4 / 16.)
    .action([](const std::string &value) { return std::stod(value); });
  program.add_argument("--tau-out")
    .help("time step size for the output")
    .default_value(1e-3)
    .action([](const std::string &value) { return std::stod(value); });
  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cout << err.what() << std::endl;
    std::cout << program;
    exit(0);
  }

  // create_for_node the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();

  mc::EmbeddedGraphReader graph_reader;
  graph_reader.append(program.get<std::string>("--input-file"), *graph);

  graph->get_vertex(0)->set_to_inflow(mc::heart_beat_inflow(program.get<double>("--heart-amplitude")));

  // set all vertices to free-outflow
  for (auto v_id : graph->get_vertex_ids())
  {
    auto v = graph->get_vertex(v_id);
    if (v->is_windkessel_outflow())
      v->set_to_free_outflow();
  }

  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, false);

  const double t_end = 1.;
  const std::size_t max_iter = 160000000;

  const auto tau = program.get<double>("--tau");
  const auto tau_out = program.get<double>("--tau-out");
  // const double tau_out = tau;
  const auto output_interval = static_cast<std::size_t>(tau_out / tau);
  std::cout << "tau = " << tau << ", tau_out = " << tau_out << ", output_interval = " << output_interval << std::endl;

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> flow_solver(MPI_COMM_WORLD, graph, dof_map_flow);
  flow_solver.set_tau(tau);
  flow_solver.use_ssp_method();

  mc::Transport transport_solver(MPI_COMM_WORLD, graph, dof_map_flow, dof_map_transport);

  std::vector<mc::Point> points;
  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<double> p_total_vertex_values;
  std::vector<double> p_static_vertex_values;

  mc::GraphFlowAndConcentrationWriter csv_writer(MPI_COMM_WORLD, "output", "data", graph, dof_map_flow, dof_map_transport);

  double total_C_e = 0;
  for (auto e_id : graph->get_active_edge_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    auto e = graph->get_edge(e_id);
    auto& data = e->get_physical_data();

    double c0 = mc::calculate_c0(data.G0, data.rho, data.A0); // [cm/s]

    // in meter
    double C_e = (1e-2*data.length) * (1e-4*data.A0) / (1060 * std::pow(1e-2*c0, 2));

    total_C_e += C_e;
  }

  // set the total flows of the leaf nodes to zero
  std::map<size_t, double> total_flows;
  for (auto v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
    if (graph->get_vertex(v_id)->is_leaf())
      total_flows[v_id] = 0.;
  }

  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    transport_solver.solve(it * tau, tau, flow_solver.get_solution());
    flow_solver.solve();

    // add total flows
    for (auto v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
      auto v = graph->get_vertex(v_id);
      if (v->is_leaf())
        total_flows[v_id] += flow_solver.get_flow_at_vessel_tip(*v) * tau;
    }

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << flow_solver.get_time() << std::endl;

      csv_writer.write(it * tau, flow_solver.get_solution(), transport_solver.get_solution());

      // double estimate parameters
      {
        // print the total flows
        double sum_of_flows = 0;
        double sum_of_out_flows = 0;
        double total_R = 1.34e8; // [Pa s / m^-3]
        double total_C = 9.45e-9; // [m^3 / Pa]
        for (auto v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
          auto v = graph->get_vertex(v_id);
          if (v->is_leaf())
          {
            std::cout << "flow at " << v_id << " = " << total_flows[v_id] << std::endl;
            sum_of_flows += total_flows[v_id];
            if (v->is_free_outflow())
            {
              sum_of_out_flows += total_flows[v_id];
            }
          }
        }
        MPI_Allreduce(&sum_of_flows, &sum_of_flows, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&sum_of_out_flows, &sum_of_out_flows, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        std::cout << "sum = " << sum_of_flows << std::endl;

        // divide resistances
        for (auto v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
          auto v = graph->get_vertex(v_id);
          if (v->is_free_outflow()) {
            auto z = total_flows[v_id] / sum_of_out_flows;
            double local_R = total_R / z;
            double R_1 = 4. / 5. * local_R;
            double R_2 = 1. / 5. * local_R;
            double C = z * (total_C - total_C_e);

            std::cout << "vertex=" << v_id << " R=" << local_R << " R_1=" << R_1 << " R_2=" << R_2 << " C=" << C << std::endl;
          }
        }
      }
    }

    // break
    if (flow_solver.get_time() > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;

  MPI_Finalize();
}
