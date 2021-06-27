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

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/set_0d_tree_boundary_conditions.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

class NonlinearLinearCoupling {
public:
  NonlinearLinearCoupling(
    std::shared_ptr<mc::GraphStorage> graph_nl,
    std::shared_ptr<mc::GraphStorage> graph_li,
    std::shared_ptr<mc::ExplicitNonlinearFlowSolver> nonlinear_solver,
    std::shared_ptr<mc::ImplicitLinearFlowSolver> linear_solver)
      : d_graph_nl(std::move(graph_nl)),
        d_graph_li(std::move(graph_li)),
        d_nonlinear_solver(std::move(nonlinear_solver)),
        d_linear_solver(std::move(linear_solver)) {}

  void add_coupled_vertices(const std::string &name) {
    add_coupled_vertices(name, name);
  }

  // TODO: calling this after having assembled the matrices will make the solvers fail
  //       -> find a way to prevent this.
  void add_coupled_vertices(const std::string &name_nl, const std::string &name_li) {
    auto v_nl = d_graph_nl->find_vertex_by_name(name_nl);
    auto v_li = d_graph_li->find_vertex_by_name(name_li);
    auto e_nl = d_graph_nl->get_edge(v_nl->get_edge_neighbors()[0]);
    auto e_li = d_graph_li->get_edge(v_li->get_edge_neighbors()[0]);
    const auto &data_nl = e_nl->get_physical_data();
    const auto &data_li = e_li->get_physical_data();
    v_nl->set_to_nonlinear_characteristic_inflow(data_li.G0, data_li.A0, data_li.rho, e_li->is_pointing_to(v_li->get_id()), 0, 0);
    v_li->set_to_linear_characteristic_inflow(mc::linear::get_C(data_nl), mc::linear::get_L(data_nl), e_nl->is_pointing_to(v_nl->get_id()), 0, 0);
    coupled_vertices.push_back({v_nl->get_id(), v_li->get_id()});
  }

  void update_linear_solver() {
    for (auto vertex_pair : coupled_vertices) {
      auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
      auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

      double p, q;
      d_nonlinear_solver->get_1d_pq_values_at_vertex(*v_nl, p, q);
      std::cout << " p=" << p << " q=" << q << std::endl;
      v_li->update_linear_characteristic_inflow(p, q);
    }
  }

  void update_nonlinear_solver() {
    for (auto vertex_pair : coupled_vertices) {
      auto v_nl = d_graph_nl->get_vertex(vertex_pair.vertex_id_1);
      auto v_li = d_graph_li->get_vertex(vertex_pair.vertex_id_2);

      double p, q;
      d_linear_solver->get_1d_pq_values_at_vertex(*v_li, p, q);
      std::cout << " p=" << p << " q=" << q << std::endl;
      v_nl->update_nonlinear_characteristic_inflow(p, q);
    }
  }

private:
  struct CoupledVertices {
    size_t vertex_id_1;
    size_t vertex_id_2;
  };

  std::vector<CoupledVertices> coupled_vertices;

private:
  std::shared_ptr<mc::GraphStorage> d_graph_nl;
  std::shared_ptr<mc::GraphStorage> d_graph_li;

  std::shared_ptr<mc::ExplicitNonlinearFlowSolver> d_nonlinear_solver;
  std::shared_ptr<mc::ImplicitLinearFlowSolver> d_linear_solver;
};

int main(int argc, char *argv[]) {
  const std::size_t degree = 1;
  const std::size_t num_micro_edges = 20;

  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const double tau = 1e-4;
    const double t_end = 0.4;
    const double tau_out = 1e-2;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);
    const size_t skip_length = 1;

    // vessel parameters
    const double vessel_length = 21.1 * 8;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double density = 1.028e-3;

    auto physical_data_1 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 9., radius, vessel_length);
    auto physical_data_2 = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length);
    // we set the viscosity to zero
    // physical_data_1.viscosity = 0;
    // physical_data_2.viscosity = 0;

    // create the ascending aorta
    auto graph_nl = std::make_shared<mc::GraphStorage>();

    auto v0_nl = graph_nl->create_vertex();
    auto v1_nl = graph_nl->create_vertex();

    auto edge_nl = graph_nl->connect(*v0_nl, *v1_nl, num_micro_edges);

    edge_nl->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1, 0, 0)}});
    edge_nl->add_physical_data(physical_data_1);

    auto graph_li = std::make_shared<mc::GraphStorage>();

    auto v0_li = graph_li->create_vertex();
    auto v1_li = graph_li->create_vertex();

    auto edge_li = graph_li->connect(*v0_li, *v1_li, num_micro_edges);

    edge_li->add_embedding_data({{mc::Point(1, 0, 0), mc::Point(2, 0, 0)}});
    edge_li->add_physical_data(physical_data_2);

    v0_nl->set_to_inflow(mc::heart_beat_inflow(4));
    v1_nl->set_name("nl_out");

    v0_li->set_name("li_in");
    // v1_li->set_to_free_outflow();
    v1_li->set_to_windkessel_outflow(1.8, 0.387);

    mc::naive_mesh_partitioner(*graph_li, PETSC_COMM_WORLD);
    mc::naive_mesh_partitioner(*graph_nl, PETSC_COMM_WORLD);

    auto dof_map_nl = std::make_shared<mc::DofMap>(graph_nl->num_vertices(), graph_nl->num_edges());
    dof_map_nl->create(PETSC_COMM_WORLD, *graph_nl, 2, degree, true);

    auto dof_map_li = std::make_shared<mc::DofMap>(graph_li->num_vertices(), graph_li->num_edges());
    dof_map_li->create(PETSC_COMM_WORLD, *graph_li, 2, degree, true);

    auto solver_nl = std::make_shared< mc::ExplicitNonlinearFlowSolver > (MPI_COMM_WORLD, graph_nl, dof_map_nl, degree);
    auto solver_li = std::make_shared< mc::ImplicitLinearFlowSolver >(PETSC_COMM_WORLD, graph_li, dof_map_li, degree);

    NonlinearLinearCoupling coupling(graph_nl, graph_li, solver_nl, solver_li);
    coupling.add_coupled_vertices("nl_out", "li_in");

    mc::GraphPVDWriter writer_li(PETSC_COMM_WORLD, "./output", "explicit_implicit_li");
    mc::GraphPVDWriter writer_nl(PETSC_COMM_WORLD, "./output", "explicit_implicit_nl");

    solver_nl->use_explicit_euler_method();
    // solver_nl.use_ssp_method();
    solver_nl->set_tau(tau);
    solver_li->setup(tau * skip_length);

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver_nl->solve();

      {
        coupling.update_linear_solver();
      }

      if (t_idx % skip_length == 0)
        solver_li->solve(tau * skip_length, t);

      {
        coupling.update_nonlinear_solver();
      }

      if (t_idx % output_interval == 0) {
        std::cout << "it = " << t_idx << std::endl;


        // linear solver
        {
          std::vector<mc::Point> points;
          std::vector<double> p_vertex_values;
          std::vector<double> q_vertex_values;
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li->p_component, solver_li->get_solution(), points, p_vertex_values);
          interpolate_to_vertices(PETSC_COMM_WORLD, *graph_li, *dof_map_li, solver_li->q_component, solver_li->get_solution(), points, q_vertex_values);

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

          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, 0, solver_nl->get_solution(), points, Q_vertex_values);
          mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, 1, solver_nl->get_solution(), points, A_vertex_values);
          mc::calculate_total_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_total_vertex_values);
          mc::calculate_static_pressure(MPI_COMM_WORLD, *graph_nl, *dof_map_nl, solver_nl->get_solution(), points, p_static_vertex_values);

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