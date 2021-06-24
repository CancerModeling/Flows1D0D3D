////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include <cmath>
#include <macrocirculation/explicit_nonlinear_flow_solver.hpp>
#include <macrocirculation/quantities_of_interest.hpp>
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

    const double tau = 1e-4;
    const double t_end = 2.;

    const size_t output_interval = 100;

    // vessel parameters
    const double vessel_length = 21.1 * 8;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double gamma = 2;
    const double density = 1.028e-3;

    auto physical_data = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, gamma, radius, vessel_length);
    // we set the viscosity to zero
    physical_data.viscosity = 0;

    // create the ascending aorta
    auto graph_nl = std::make_shared<mc::GraphStorage>();

    auto v0_nl = graph_nl ->create_vertex();
    auto v1_nl = graph_nl ->create_vertex();

    auto edge_nl = graph_nl->connect(*v0_nl, *v1_nl, num_micro_edges);

    edge_nl->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1, 0, 0)}});
    edge_nl->add_physical_data(physical_data);

    auto graph_li = std::make_shared<mc::GraphStorage>();

    auto v0_li = graph_li ->create_vertex();
    auto v1_li = graph_li ->create_vertex();

    auto edge_li = graph_li->connect(*v0_li, *v1_li, num_micro_edges);

    edge_li->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(2, 0, 0)}});
    edge_li->add_physical_data(physical_data);

    v0_nl->set_to_inflow(mc::heart_beat_inflow(0.4));
    v1_nl->set_to_nonlinear_characteristic_inflow(physical_data.G0, physical_data.A0, physical_data.rho, false, 0, 0);

    v0_li->set_to_linear_characteristic_inflow(mc::LinearFlowSolver::get_C(*edge_li), mc::LinearFlowSolver::get_L(*edge_li), true, 0, 0);
    v1_li->set_to_free_outflow();

    mc::naive_mesh_partitioner(*graph_li, PETSC_COMM_WORLD);
    mc::naive_mesh_partitioner(*graph_nl, PETSC_COMM_WORLD);

    auto dof_map_nl = std::make_shared<mc::DofMap>(graph_nl->num_vertices(), graph_nl->num_edges());
    dof_map_nl->create(PETSC_COMM_WORLD, *graph_nl, 2, degree, true);

    auto dof_map_li = std::make_shared<mc::DofMap>(graph_li->num_vertices(), graph_li->num_edges());
    dof_map_li->create(PETSC_COMM_WORLD, *graph_li, 2, degree, true);

    mc::ExplicitNonlinearFlowSolver<degree> solver_nl(MPI_COMM_WORLD, graph_nl, dof_map_nl);
    solver_nl.use_explicit_euler_method();
    solver_nl.set_tau(tau);

    mc::LinearFlowSolver solver_li(PETSC_COMM_WORLD, graph_li, dof_map_li, degree);
    solver_li.setup(tau);

    mc::GraphPVDWriter writer_li(PETSC_COMM_WORLD, "./output", "explicit_implicit_li");
    mc::GraphPVDWriter writer_nl(PETSC_COMM_WORLD, "./output", "explicit_implicit_nl");

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver_nl.solve();

      {
        double Q,A;
        solver_nl.get_1d_values_at_vertex(*v1_nl, Q, A);
        const auto& param = v0_li->get_linear_characteristic_data();
        const double p = physical_data.G0 * (std::sqrt(A/physical_data.A0) - 1);
        v0_li->update_linear_characteristic_inflow(p, Q);

        std::cout << " p=" << p << " Q=" << Q << std::endl;
      }

      solver_li.solve(tau, t);

      {
        double p,q;
        solver_li.get_1d_values_at_vertex(*v0_li, p, q);
        const auto& param = v1_nl->get_nonlinear_characteristic_data();
        v1_nl->update_nonlinear_characteristic_inflow(p, q);

        std::cout << " p=" << p << " q=" << q << std::endl;
      }

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;


        // linear solver
        {
          std::vector<mc::Point> points;
          std::vector<double> p_vertex_values;
          std::vector<double> q_vertex_values;
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li.p_component, solver_li.get_solution(), points, p_vertex_values);
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li.q_component, solver_li.get_solution(), points, q_vertex_values);

          writer_li.set_points(points);
          writer_li.add_vertex_data("p", p_vertex_values);
          writer_li.add_vertex_data("q", q_vertex_values);
          writer_li.write(t);
        }

        // nonlinear solver
        {
          std::vector<mc::Point> points;
          std::vector<double> Q_vertex_values;
          std::vector<double> A_vertex_values;
          std::vector<double> p_total_vertex_values;
          std::vector<double> p_static_vertex_values;

          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, 0, solver_nl.get_solution(), points, Q_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, 1, solver_nl.get_solution(), points, A_vertex_values);
          mc::calculate_total_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl.get_solution(), points, p_total_vertex_values);
          mc::calculate_static_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl.get_solution(), points, p_static_vertex_values);

          writer_nl.set_points(points);
          writer_nl.add_vertex_data("Q", Q_vertex_values);
          writer_nl.add_vertex_data("A", A_vertex_values);
          writer_nl.add_vertex_data("p_static", p_static_vertex_values);
          writer_nl.add_vertex_data("p_total", p_total_vertex_values);
          writer_nl.write(t);
        }
      }
    }
  }

  CHKERRQ(PetscFinalize());
}