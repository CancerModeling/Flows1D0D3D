////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <cmath>
#include <memory>
#include <petsc.h>
#include <utility>
#include <cxxopts.hpp>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/set_0d_tree_boundary_conditions.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;


int main(int argc, char *argv[]) {
  const std::size_t degree = 2;

  cxxopts::Options options(argv[0], "Linearized solver for breast geometry");
  options.add_options()                                                                                                             //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output/"))                     //
    ("tau", "time step size", cxxopts::value<double>()->default_value("1e-4"))                                //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                   //
    ("t-end", "Endtime for simulation", cxxopts::value<double>()->default_value("10"))                                               //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if (!args.unmatched().empty()) {
    std::cout << "The following arguments were unmatched: " << std::endl;
    for (auto &it : args.unmatched())
      std::cout << " " << it;
    std::cout << "\nAre they part of petsc or a different auxillary library?" << std::endl;
  }

  const double t_end = args["t-end"].as<double>();
  const std::size_t max_iter = 160000000;

  const auto tau = args["tau"].as<double>();
  const auto tau_out = args["tau-out"].as<double>();

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    // vessel parameters

    // create the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();

    mc::EmbeddedGraphReader graph_reader;
    graph_reader.append("./data/meshes/coarse-network-geometry.json", *graph);

    graph->find_vertex_by_name("bg_132")->set_to_inflow(mc::heart_beat_inflow(0.035));
    graph->find_vertex_by_name("bg_135")->set_to_inflow(mc::heart_beat_inflow(0.035));
    graph->find_vertex_by_name("bg_141")->set_to_inflow(mc::heart_beat_inflow(0.035));
    graph->find_vertex_by_name("bg_119")->set_to_inflow(mc::heart_beat_inflow(0.035));

    mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

    mc::set_0d_tree_boundary_conditions(graph, "bg_");

    auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map->create(PETSC_COMM_WORLD, *graph, 2, degree, true);

    mc::LinearFlowSolver solver(PETSC_COMM_WORLD, graph, dof_map, degree);
    solver.setup(tau);

    mc::GraphPVDWriter writer(MPI_COMM_WORLD, args["output-directory"].as<std::string>(), "breast_geometry_linearized");

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver.solve(tau, t);

      if (t_idx % output_interval == 0) {
        std::cout << "iter = " << t_idx << ", t = " << t << std::endl;

        std::vector<mc::Point> points;
        std::vector<double> p_vertex_values;
        std::vector<double> q_vertex_values;
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.p_component, solver.get_solution(), points, p_vertex_values);
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.q_component, solver.get_solution(), points, q_vertex_values);

        writer.set_points(points);
        writer.add_vertex_data("p", p_vertex_values);
        writer.add_vertex_data("q", q_vertex_values);
        writer.write(t);

        for (auto &v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
          auto &vertex = *graph->get_vertex(v_id);
          if (vertex.is_vessel_tree_outflow()) {
            const auto &vertex_dof_map = dof_map->get_local_dof_map(vertex);
            const auto &vertex_dofs = vertex_dof_map.dof_indices();
            const auto &u = solver.get_solution();

            std::cout << "rank = " << mc::mpi::rank(MPI_COMM_WORLD) << std::endl;
            std::cout << "vertex = " << vertex.get_name() << "\n";
            std::cout << " p = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << u.get(vertex_dofs[k]) << ", ";
            std::cout << std::endl;
            std::cout << " R = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << vertex.get_vessel_tree_data().resistances[k] << ", ";
            std::cout << std::endl;
            std::cout << " C = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << vertex.get_vessel_tree_data().capacitances[k] << ", ";
            std::cout << std::endl;
          } else {
            // std::cout << "vertex = " << vertex.get_name() << " (no tree outflow)\n";
          }
        }

        std::cout << "iter = " << t_idx << ", t = " << t << std::endl;
      }
    }
  }

  CHKERRQ(PetscFinalize());
}