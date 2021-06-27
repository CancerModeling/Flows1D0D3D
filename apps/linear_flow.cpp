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
  const std::size_t num_micro_edges = 20;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const double tau = 1e-3;
    const double t_end = 1;

    const size_t output_interval = 100;

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
    auto edge1 = graph->connect(*v0, *v1, num_micro_edges);

    auto physical_data = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
    physical_data.viscosity = 0.;

    edge1->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(2, 0, 0)}});
    edge1->add_physical_data(physical_data);

    const double p_in = 5.;
    const double q_in = 4.;

    v0->set_to_linear_characteristic_inflow(mc::LinearFlowSolver::get_C(*edge1), mc::LinearFlowSolver::get_L(*edge1),true, p_in, q_in);
    v1->set_to_linear_characteristic_inflow(mc::LinearFlowSolver::get_C(*edge1), mc::LinearFlowSolver::get_L(*edge1),false, p_in, q_in);

    mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);

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

        {
          double p,q;
          solver.evaluate_1d_pq_values(*edge1, 0, p, q);
          std::cout << std::abs(p - p_in) << " " << std::abs(q - q_in) << " ";
          solver.evaluate_1d_pq_values(*edge1, 0.5, p, q);
          std::cout << std::abs(p - p_in) << " " << std::abs(q - q_in) << " ";
          solver.evaluate_1d_pq_values(*edge1, 1, p, q);
          std::cout << std::abs(p - p_in) << " " << std::abs(q - q_in) << " ";
          std::cout << std::endl;
        }

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