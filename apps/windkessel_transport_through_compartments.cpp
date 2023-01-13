////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <memory>
#include <utility>

#include "petsc.h"

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/csv_vessel_tip_writer.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/linearized_flow_upwind_evaluator.hpp"
#include "macrocirculation/petsc/petsc_mat.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  cxxopts::Options options(argv[0], "Implicit transport solver with artificial flow.");
  options.add_options()                                                                                                           //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(1e-3)))                     //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                 //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("6"))                                   //
    ("degree", "The FEM degree of the transport", cxxopts::value<size_t>()->default_value("2"))                                   //
    ("num-micro-edges", "The number of micro edges on each vessel (we have two)", cxxopts::value<size_t>()->default_value("40"))  //
    ("backward-orientation", "Should the orientation of the vessels be flipped?", cxxopts::value<bool>()->default_value("false")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves implicit transport problem with artificial flow."));
  {
    const double tau = args["tau"].as<double>();
    const double tau_out = args["tau-out"].as<double>();
    const double t_end = args["t-end"].as<double>();
    const size_t degree = args["degree"].as<size_t>();
    const size_t num_micro_edges = args["num-micro-edges"].as<size_t>();

    const bool backward_orientation = args["backward-orientation"].as<bool>();

    // vessel parameters
    const double vessel_length = 1.;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double density = 1.028e-3;

    auto physical_data_short = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
    auto physical_data_long = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);


    // create_for_node the geometry of the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();
    auto v0 = graph->create_vertex();
    auto v1 = graph->create_vertex();
    auto v2 = graph->create_vertex();

    std::shared_ptr<mc::Edge> edge_1;
    std::shared_ptr<mc::Edge> edge_2;
    if (backward_orientation) {
      edge_1 = graph->connect(*v1, *v0, num_micro_edges);
      edge_1->add_embedding_data({{mc::Point(0.5, 0, 0), mc::Point(0, 0, 0)}});
      edge_2 = graph->connect(*v2, *v1, num_micro_edges);
      edge_2->add_embedding_data({{mc::Point(1., 0, 0), mc::Point(0.5, 0, 0)}});
    } else {
      edge_1 = graph->connect(*v0, *v1, num_micro_edges);
      edge_1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(0.5, 0, 0)}});
      edge_2 = graph->connect(*v1, *v2, num_micro_edges);
      edge_2->add_embedding_data({{mc::Point(0.5, 0, 0), mc::Point(1., 0, 0)}});
    }
    edge_2->add_physical_data(physical_data_long);
    edge_1->add_physical_data(physical_data_short);

    // v2->set_to_vessel_tree_outflow(5 * (133.333) * 1e-2, {1., 1., 1.}, {1., 1., 1.}, 1);
    v2->set_to_vessel_tree_outflow(5 * (133.333) * 1e-2, {1., 1.}, {1., 1.}, {radius, radius}, 3);
    v0->set_to_inflow_with_fixed_flow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });

    graph->finalize_bcs();

    // partition graph
    mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

    const std::size_t max_iter = 1600000;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    // configure solver
    auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    mc::DofMap::create(MPI_COMM_WORLD, {graph}, {dof_map_transport}, 1, degree, [](auto, const mc::Vertex &v) -> size_t {
      if (v.is_windkessel_outflow())
        return 1;
      else if (v.is_vessel_tree_outflow())
        return v.get_vessel_tree_data().resistances.size();
      return 0;
    });

    auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, true);

    auto upwind_evaluator = std::make_shared<mc::LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
    auto flow_inside = [backward_orientation](double t, const mc::Point &/*p*/) -> double {
      const double sgn = backward_orientation ? -1. : +1.;
      if (t < 2)
        return sgn;
      else
        return 0;
    };
    auto flow_outside = [&dof_map_transport](double t, const mc::Vertex &v, std::vector<double> &p_c) {
      auto num_dof = dof_map_transport->get_local_dof_map(v).num_local_dof();
      p_c.resize(num_dof + 1);
      for (size_t k = 0; k < num_dof + 1; k += 1)
        p_c[k] = 0;
      // p_c[p_c.size()-2] = p_c[p_c.size()-1];
      if ( t > 2)
      {
        p_c[0] = 2;
        p_c[1] = 1;
        p_c[2] = 1;
      }
      if ( t > 4)
      {
        p_c[0] = 1;
        p_c[1] = 2;
        p_c[2] = 2;
      }
    };
    auto variable_upwind_provider = std::make_shared<mc::EmbeddedUpwindProvider>(graph, flow_inside, flow_outside);

    mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

    transport_solver.set_inflow_function([](double t) -> double { return (t < 0.25) ? 1 : 0; });

    mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

    mc::CSVVesselTipWriter vessel_tip_writer(MPI_COMM_WORLD,
                                             "output", "implicit_transport_with_0d",
                                             graph,
                                             {dof_map_transport, transport_solver.get_dof_maps_volume().front()},
                                             {"c", "V"});

    std::vector<mc::Point> points;

    const auto begin_t = std::chrono::steady_clock::now();
    double t = 0;
    for (std::size_t it = 0; it < max_iter; it += 1) {
      transport_solver.solve(tau, t + tau);
      t += tau;

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", time = " << t << std::endl;

        // save solution
        std::vector<double> c_vertex_values;
        std::vector<double> v_vertex_values;
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *variable_upwind_provider, t, points, v_vertex_values);

        pvd_writer.set_points(points);
        pvd_writer.add_vertex_data("c", c_vertex_values);
        pvd_writer.add_vertex_data("v", v_vertex_values);
        pvd_writer.write(t);

        vessel_tip_writer.write(t, {transport_solver.get_solution(), transport_solver.get_volumes()});
      }

      // break
      if (t > t_end + 1e-12)
        break;
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}
