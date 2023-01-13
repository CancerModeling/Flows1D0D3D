////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  const std::size_t degree = 2;
  const std::size_t num_micro_edges = 10;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(MPI_COMM_WORLD) << std::endl;

    const double tau = 2.5e-4 / 32.;
    // const double t_end = 2.5e-4 / 16 * 100;
    const double t_end = 4;
    const double tau_out = 1e-2;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    // vessel parameters
    const double vessel_length = 2.5;
    // const double radius = 0.403;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double density = 1.028e-3;

    auto physical_data_1 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 9., radius, vessel_length / 2.);
    auto physical_data_2 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
    // we set the viscosity to zero
    // physical_data_1.viscosity = 0;
    // physical_data_2.viscosity = 0;

    bool nl_forward = false;
    bool li_forward = false;

    // create the ascending aorta
    auto graph_nl = std::make_shared<mc::GraphStorage>();

    auto v0_nl = graph_nl->create_vertex();
    auto v1_nl = graph_nl->create_vertex();
    auto v2_nl = graph_nl->create_vertex();


    auto edge_0_nl = graph_nl->connect(*v0_nl, *v1_nl, num_micro_edges);
    std::shared_ptr<mc::Edge> edge_1_nl;
    if (nl_forward)
      edge_1_nl = graph_nl->connect(*v1_nl, *v2_nl, num_micro_edges);
    else
      edge_1_nl = graph_nl->connect(*v2_nl, *v1_nl, num_micro_edges);

    edge_0_nl->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(0.5, 0, 0)}});
    edge_0_nl->add_physical_data(physical_data_1);
    if (nl_forward)
      edge_1_nl->add_embedding_data({{mc::Point(0.5, 0, 0), mc::Point(1, 0, 0)}});
    else
      edge_1_nl->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(0.5, 0, 0)}});
    edge_1_nl->add_physical_data(physical_data_1);

    auto graph_li = std::make_shared<mc::GraphStorage>();

    auto v0_li = graph_li->create_vertex();
    auto v1_li = graph_li->create_vertex();
    auto v2_li = graph_li->create_vertex();

    std::shared_ptr<mc::Edge> edge_0_li;
    if (li_forward)
      edge_0_li = graph_li->connect(*v0_li, *v1_li, num_micro_edges);
    else
      edge_0_li = graph_li->connect(*v1_li, *v0_li, num_micro_edges);
    auto edge_1_li = graph_li->connect(*v2_li, *v1_li, num_micro_edges);

    if (li_forward)
      edge_0_li->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(1.5, 0, 0)}});
    else
      edge_0_li->add_embedding_data({{mc::Point(1.5, 0, 0), mc::Point(1, 0, 0)}});
    edge_0_li->add_physical_data(physical_data_2);
    edge_1_li->add_embedding_data({{mc::Point(2., 0, 0), mc::Point(1.5, 0, 0)}});
    edge_1_li->add_physical_data(physical_data_2);

    v0_nl->set_to_inflow_with_fixed_flow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });
    v2_nl->set_name("nl_out");

    v0_li->set_name("li_in");
    // v1_li->set_to_free_outflow();

    // v2_li->set_to_free_outflow();
    v2_li->set_to_vessel_tree_outflow(5.0 * 1.333322, {1.8}, {3870}, {1.}, 1.);

    // v2_li->set_name("windkessel_outflow");
    // mc::set_0d_tree_boundary_conditions(graph_li, "windkessel_outflow");

    mc::naive_mesh_partitioner(*graph_li, MPI_COMM_WORLD);
    mc::naive_mesh_partitioner(*graph_nl, MPI_COMM_WORLD);

    auto coupling = std::make_shared<mc::NonlinearLinearCoupling>(MPI_COMM_WORLD, graph_nl, graph_li);
    coupling->add_coupled_vertices("nl_out", "li_in");

    mc::CoupledExplicitImplicit1DSolver solver(MPI_COMM_WORLD, coupling, graph_nl, graph_li, degree, degree);
    solver.get_explicit_solver()->use_explicit_euler_method();

    auto dof_map_nl = solver.get_explicit_dof_map();
    auto dof_map_li = solver.get_implicit_dof_map();

    auto dof_map_transport_nl = std::make_shared<mc::DofMap>(*graph_nl);
    auto dof_map_transport_li = std::make_shared<mc::DofMap>(*graph_li);
    mc::DofMap::create(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, 1, degree, [](const mc::GraphStorage &, const mc::Vertex &v) {
      if (v.is_vessel_tree_outflow())
        return 1;
      return 0;
    });

    std::cout << mc::mpi::rank(MPI_COMM_WORLD) << " last global dof " << dof_map_transport_nl->last_global_dof() << std::endl;
    for (auto eid : graph_nl->get_edge_ids()) {
      auto edge = graph_nl->get_edge(eid);
      auto ldof_map = dof_map_transport_nl->get_local_dof_map(*edge);
      std::vector<size_t> dof_indices(ldof_map.num_basis_functions());
      for (size_t meid = 0; meid < edge->num_micro_edges(); meid += 1) {
        ldof_map.dof_indices(meid, 0, dof_indices);
        std::cout << "[" << mc::mpi::rank(MPI_COMM_WORLD) << "] nonlinear "
                  << "macro-edge = " << eid << ", micro-edge = " << meid << ", dof-indices " << dof_indices << std::endl;
      }
    }
    for (auto vid : graph_nl->get_vertex_ids()) {
      auto vertex = graph_nl->get_vertex(vid);
      if (!vertex->is_leaf())
        continue;
      auto ldof_map = dof_map_transport_nl->get_local_dof_map(*vertex);
      std::cout << "[" << mc::mpi::rank(MPI_COMM_WORLD) << "] nonlinear "
                << "macro-vertex= " << vid << ", dof-indices " << ldof_map.dof_indices() << std::endl;
    }

    std::cout << mc::mpi::rank(MPI_COMM_WORLD) << " last global dof " << dof_map_transport_li->last_global_dof() << std::endl;
    for (auto eid : graph_li->get_edge_ids()) {
      auto edge = graph_li->get_edge(eid);
      auto ldof_map = dof_map_transport_li->get_local_dof_map(*edge);
      std::vector<size_t> dof_indices(ldof_map.num_basis_functions());
      for (size_t meid = 0; meid < edge->num_micro_edges(); meid += 1) {
        ldof_map.dof_indices(meid, 0, dof_indices);
        std::cout << "[" << mc::mpi::rank(MPI_COMM_WORLD) << "] linear "
                  << "macro-edge = " << eid << ", micro-edge = " << meid << ", dof-indices " << dof_indices << std::endl;
      }
    }
    for (auto vid : graph_li->get_vertex_ids()) {
      auto vertex = graph_li->get_vertex(vid);
      if (!vertex->is_leaf())
        continue;
      auto ldof_map = dof_map_transport_li->get_local_dof_map(*vertex);
      std::cout << "[" << mc::mpi::rank(MPI_COMM_WORLD) << "] linear "
                << "macro-vertex= " << vid << ", dof-indices " << ldof_map.dof_indices() << std::endl;
    }

    auto upwind_evaluator_nl = std::make_shared<mc::NonlinearFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_nl, dof_map_nl);
    auto variable_upwind_provider_nl = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator_nl, solver.get_explicit_solver());

    auto upwind_evaluator_li = std::make_shared<mc::LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph_li, dof_map_li);
    auto variable_upwind_provider_li = std::make_shared<mc::UpwindProviderLinearizedFlow>(graph_li, upwind_evaluator_li, solver.get_implicit_solver());

    mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, {graph_nl, graph_li}, {dof_map_transport_nl, dof_map_transport_li}, {variable_upwind_provider_nl, variable_upwind_provider_li}, degree);

    transport_solver.set_inflow_function([](double t) { return 1.; });

    auto solver_nl = solver.get_explicit_solver();
    auto solver_li = solver.get_implicit_solver();

    mc::GraphPVDWriter writer_li(MPI_COMM_WORLD, "./output", "explicit_implicit_li");
    mc::GraphPVDWriter writer_nl(MPI_COMM_WORLD, "./output", "explicit_implicit_nl");

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

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;

        std::cout << "norm transport solution: " << transport_solver.get_solution().norm2() << std::endl;
        std::cout << "norm transport rhs: " << transport_solver.get_rhs().norm2() << std::endl;

        // transport_solver.get_solution().print();
        // transport_solver.get_rhs().print();
        // transport_solver.get_mat().print();

        // linear solver
        {
          std::vector<double> p_vertex_values;
          std::vector<double> q_vertex_values;
          std::vector<double> c_vertex_values;
          std::vector<double> v_vertex_values;
          std::vector<double> A_vertex_values;
          interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, solver_li->p_component, solver_li->get_solution(), points, p_vertex_values);
          interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_li, solver_li->q_component, solver_li->get_solution(), points, q_vertex_values);
          interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *dof_map_transport_li, 0, transport_solver.get_solution(), points, c_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_li, *variable_upwind_provider_li, t, points, v_vertex_values);
          auto trafo = [](double p, const mc::Edge &e) {
            return e.get_physical_data().A0 + mc::linear::get_C(e.get_physical_data()) * p;
          };
          mc::interpolate_transformation(MPI_COMM_WORLD, *graph_li, *dof_map_li, solver_li->p_component, solver_li->get_solution(), trafo, points, A_vertex_values);

          writer_li.set_points(points);
          writer_li.add_vertex_data("p", p_vertex_values);
          writer_li.add_vertex_data("q", q_vertex_values);
          writer_li.add_vertex_data("c", c_vertex_values);
          writer_li.add_vertex_data("A", A_vertex_values);
          writer_li.add_vertex_data("v", v_vertex_values);
          writer_li.write(t);
        }

        // nonlinear solver
        {
          std::vector<double> Q_vertex_values;
          std::vector<double> A_vertex_values;
          std::vector<double> p_total_vertex_values;
          std::vector<double> p_static_vertex_values;
          std::vector<double> c_vertex_values;
          std::vector<double> v_vertex_values;

          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->Q_component, solver_nl->get_solution(), points, Q_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->A_component, solver_nl->get_solution(), points, A_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_transport_nl, 0, transport_solver.get_solution(), points, c_vertex_values);
          mc::calculate_total_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_total_vertex_values);
          mc::calculate_static_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_static_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *variable_upwind_provider_nl, t, points, v_vertex_values);

          // VecView(transport_solver.get_solution().get_vec(), PETSC_VIEWER_STDOUT_WORLD);
          writer_nl.set_points(points);
          writer_nl.add_vertex_data("Q", Q_vertex_values);
          writer_nl.add_vertex_data("A", A_vertex_values);
          writer_nl.add_vertex_data("p_static", p_static_vertex_values);
          writer_nl.add_vertex_data("p_total", p_total_vertex_values);
          writer_nl.add_vertex_data("c", c_vertex_values);
          writer_nl.add_vertex_data("v", v_vertex_values);
          writer_nl.write(t);
        }
      }
    }
  }

  CHKERRQ(PetscFinalize());
}