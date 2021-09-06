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

#include "macrocirculation/0d_boundary_conditions.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/coupled_explicit_implicit_1d_solver.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/linearized_flow_upwind_evaluator.hpp"
#include "macrocirculation/nonlinear_flow_upwind_evaluator.hpp"
#include "macrocirculation/nonlinear_linear_coupling.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 10;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const double tau = 2.5e-4 / 32.;
    const double t_end = 3.;
    const double tau_out = 1e-2;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);
    const size_t skip_length = 1;

    // vessel parameters
    const double vessel_length = 2.5;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double density = 1.028e-3;

    auto physical_data_1 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 9., radius, vessel_length / 2.);
    auto physical_data_2 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
    // we set the viscosity to zero
    // physical_data_1.viscosity = 0;
    // physical_data_2.viscosity = 0;

    // create the ascending aorta
    auto graph_nl = std::make_shared<mc::GraphStorage>();

    auto v0_nl = graph_nl->create_vertex();
    auto v1_nl = graph_nl->create_vertex();
    auto v2_nl = graph_nl->create_vertex();

    auto edge_0_nl = graph_nl->connect(*v0_nl, *v1_nl, num_micro_edges);
    // auto edge_1_nl = graph_nl->connect(*v1_nl, *v2_nl, num_micro_edges);
    auto edge_1_nl = graph_nl->connect(*v2_nl, *v1_nl, num_micro_edges);

    edge_0_nl->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(0.5, 0, 0)}});
    edge_0_nl->add_physical_data(physical_data_1);
    // edge_1_nl->add_embedding_data({{mc::Point(0.5, 0, 0), mc::Point(1, 0, 0)}});
    edge_1_nl->add_embedding_data({{mc::Point(1.0, 0, 0), mc::Point(0.5, 0, 0)}});
    edge_1_nl->add_physical_data(physical_data_1);

    auto graph_li = std::make_shared<mc::GraphStorage>();

    auto v0_li = graph_li->create_vertex();
    auto v1_li = graph_li->create_vertex();
    auto v2_li = graph_li->create_vertex();

    auto edge_0_li = graph_li->connect(*v0_li, *v1_li, num_micro_edges);
    // auto edge_1_li = graph_li->connect(*v1_li, *v2_li, num_micro_edges);
    auto edge_1_li = graph_li->connect(*v2_li, *v1_li, num_micro_edges);

    edge_0_li->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1.5, 0, 0)}});
    edge_0_li->add_physical_data(physical_data_2);
    //edge_1_li->add_embedding_data({{mc::Point(1.5, 0, 0), mc::Point(2, 0, 0)}});
    edge_1_li->add_embedding_data({{mc::Point(2, 0, 0), mc::Point(1.5, 0, 0)}});
    edge_1_li->add_physical_data(physical_data_2);

    v0_nl->set_to_inflow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });
    v2_nl->set_name("nl_out");

    v0_li->set_name("li_in");
    // v1_li->set_to_free_outflow();

    //v2_li->set_to_windkessel_outflow(1.8, 0.387);
    // v2_li->set_to_free_outflow();
    v2_li->set_to_windkessel_outflow(18, 387);

    // v2_li->set_name("windkessel_outflow");
    // mc::set_0d_tree_boundary_conditions(graph_li, "windkessel_outflow");

    mc::naive_mesh_partitioner(*graph_li, PETSC_COMM_WORLD);
    mc::naive_mesh_partitioner(*graph_nl, PETSC_COMM_WORLD);

    auto coupling = std::make_shared<mc::NonlinearLinearCoupling>(MPI_COMM_WORLD, graph_nl, graph_li);
    coupling->add_coupled_vertices("nl_out", "li_in");

    mc::CoupledExplicitImplicit1DSolver solver(MPI_COMM_WORLD, coupling, graph_nl, graph_li, degree, degree);
    solver.get_explicit_solver()->use_ssp_method();

    auto dof_map_nl = solver.get_explicit_dof_map();
    auto dof_map_li = solver.get_implicit_dof_map();

    auto dof_map_transport_nl = std::make_shared<mc::DofMap>(*graph_nl);
    dof_map_transport_nl->create(PETSC_COMM_WORLD, *graph_nl, 1, degree, 0, true);
    auto dof_map_transport_li = std::make_shared<mc::DofMap>(*graph_li);
    dof_map_transport_li->create(PETSC_COMM_WORLD, *graph_li, 1, degree, dof_map_transport_nl->last_global_dof() + 1, true);

    auto upwind_evaluator_nl = std::make_shared<mc::NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_nl, dof_map_nl);
    auto variable_upwind_provider_nl = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator_nl, solver.get_explicit_solver());

    auto upwind_evaluator_li = std::make_shared<mc::LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_li, dof_map_li);
    auto variable_upwind_provider_li = std::make_shared<mc::UpwindProviderLinearizedFlow>(graph_li, upwind_evaluator_li, solver.get_implicit_solver());

    mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, {variable_upwind_provider_nl, variable_upwind_provider_li}, degree);

    auto solver_nl = solver.get_explicit_solver();
    auto solver_li = solver.get_implicit_solver();

    mc::GraphPVDWriter writer_li(PETSC_COMM_WORLD, "./output", "explicit_implicit_li");
    mc::GraphPVDWriter writer_nl(PETSC_COMM_WORLD, "./output", "explicit_implicit_nl");

    std::vector<mc::Point> points;
    std::vector<double> vessel_A0_li;
    mc::fill_with_vessel_A0(MPI_COMM_WORLD, *graph_li, points, vessel_A0_li);

    solver.setup(tau);

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      solver.solve(tau, t);
      upwind_evaluator_nl->init(t + tau, solver.get_explicit_solver()->get_solution());
      upwind_evaluator_li->init(t + tau, solver.get_implicit_solver()->get_solution());
      transport_solver.solve(tau, t + tau);

      t += tau;

      if (t_idx == 95478)
      {
        std::cout << "it = " << t_idx << " t  = " << t << std::endl;
        std::cout << transport_solver.get_solution() << std::endl;
        std::cout << transport_solver.get_mat() << std::endl;
        std::cout << transport_solver.get_rhs() << std::endl;
      }
      if (transport_solver.get_solution().norm2() < 1e-10 && t > 0.1)
      {
        std::cout << "it = " << t_idx << " t  = " << t << std::endl;
        std::cout << transport_solver.get_solution() << std::endl;
        std::cout << transport_solver.get_mat() << std::endl;
        std::cout << transport_solver.get_rhs() << std::endl;
        return 0;
      }

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;

        // linear solver
        {
          std::vector<double> p_vertex_values;
          std::vector<double> q_vertex_values;
          std::vector<double> c_vertex_values;
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li->p_component, solver_li->get_solution(), points, p_vertex_values);
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li->q_component, solver_li->get_solution(), points, q_vertex_values);
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_transport_li, 0, transport_solver.get_solution(), points, c_vertex_values);


          writer_li.set_points(points);
          writer_li.add_vertex_data("p", p_vertex_values);
          writer_li.add_vertex_data("q", q_vertex_values);
          writer_li.add_vertex_data("c", c_vertex_values);
          writer_li.add_vertex_data("A", vessel_A0_li);
          writer_li.write(t);
        }

        // nonlinear solver
        {
          std::vector<double> Q_vertex_values;
          std::vector<double> A_vertex_values;
          std::vector<double> p_total_vertex_values;
          std::vector<double> p_static_vertex_values;
          std::vector<double> c_vertex_values;

          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->Q_component, solver_nl->get_solution(), points, Q_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->A_component, solver_nl->get_solution(), points, A_vertex_values);
          mc::interpolate_to_vertices(PETSC_COMM_WORLD, *graph_nl, *dof_map_transport_nl, 0, transport_solver.get_solution(), points, c_vertex_values);
          mc::calculate_total_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_total_vertex_values);
          mc::calculate_static_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_static_vertex_values);

          writer_nl.set_points(points);
          writer_nl.add_vertex_data("Q", Q_vertex_values);
          writer_nl.add_vertex_data("A", A_vertex_values);
          writer_nl.add_vertex_data("p_static", p_static_vertex_values);
          writer_nl.add_vertex_data("p_total", p_total_vertex_values);
          writer_nl.add_vertex_data("c", c_vertex_values);
          writer_nl.write(t);
        }
      }
    }
  }

  CHKERRQ(PetscFinalize());
}