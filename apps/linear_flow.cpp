////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include <cmath>
#include <macrocirculation/set_0d_tree_boundary_conditions.hpp>
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;


int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 44;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const double tau = 1e-3;
    const double t_end = 1.;

    const size_t output_interval = 1;

    // vessel parameters
    const double vessel_length = 42.2;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double gamma = 9;
    const double density = 1.028e-3;

    // create the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();

    auto v0 = graph->create_vertex();
    auto v1 = graph->create_vertex();
    auto v2 = graph->create_vertex();
    auto v3 = graph->create_vertex();
    auto v4 = graph->create_vertex();
    auto edge1 = graph->connect(*v0, *v1, num_micro_edges);
    auto edge2 = graph->connect(*v1, *v2, num_micro_edges);
    auto edge3 = graph->connect(*v1, *v3, num_micro_edges);
    auto edge4 = graph->connect(*v1, *v4, num_micro_edges);

    v0->set_to_inflow(mc::heart_beat_inflow(4.));
    //v2->set_to_free_outflow();
    // v1->set_to_vessel_tree_outflow(5.0 * 1.333322, {1.8, 1.8, 1.8}, {0.387, 0.387, 0.387});
    // v1->set_to_windkessel_outflow(1.8, 0.387);
    v2->set_to_windkessel_outflow(1.8, 0.387);
    v3->set_to_windkessel_outflow(1.8, 0.387);
    v4->set_to_windkessel_outflow(1.8, 0.387);

    v2->set_name("b_2");
    v3->set_name("b_3");
    v4->set_name("b_4");

    auto physical_data_1 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
    auto physical_data_2 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
    auto physical_data_3 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
    auto physical_data_4 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);

    edge1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1, 0, 0)}});
    edge1->add_physical_data(physical_data_1);
    edge2->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(2, 0, 0)}});
    edge2->add_physical_data(physical_data_2);
    edge3->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1, -1, 0)}});
    edge3->add_physical_data(physical_data_3);
    edge4->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1, +1, 0)}});
    edge4->add_physical_data(physical_data_4);

    mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

    // mc::set_0d_tree_boundary_conditions(graph, "b_");

    auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map->create(PETSC_COMM_WORLD, *graph, 2, degree, true);

    mc::LinearFlowSolver solver(PETSC_COMM_WORLD, graph, dof_map, degree);
    solver.setup(tau);

    mc::GraphPVDWriter writer(PETSC_COMM_WORLD, "./output", "linear_flow");

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver.solve(tau, t);

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;

        std::vector<mc::Point> points;
        std::vector<double> p_vertex_values;
        std::vector<double> q_vertex_values;
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.p_component, solver.get_solution(), points, p_vertex_values);
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.q_component, solver.get_solution(), points, q_vertex_values);

        writer.set_points(points);
        writer.add_vertex_data("p", p_vertex_values);
        writer.add_vertex_data("q", q_vertex_values);
        writer.write(t);
      }
    }
  }

  CHKERRQ(PetscFinalize());
}