////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <memory>
#include <utility>

#include "petsc.h"

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/transport.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/right_hand_side_evaluator.hpp"


namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

class ConstantUpwindProvider : public mc::UpwindProvider {
public:
  ConstantUpwindProvider(double speed)
      : d_speed(speed) {}

  ~ConstantUpwindProvider() override = default;

  void get_values_at_qp(double t,
                        const mc::Edge &edge,
                        size_t micro_edge,
                        const mc::QuadratureFormula &qf,
                        std::vector<double> &v_qp) override {
    assert(qf.size() == v_qp.size());

    for (size_t k = 0; k < qf.size(); k += 1)
      v_qp[k] = d_speed;
  }

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const mc::Edge &edge,
                           std::vector<double> &v_qp) override {
    assert(v_qp.size() == edge.num_micro_vertices());

    for (size_t k = 0; k < edge.num_micro_vertices(); k += 1)
      v_qp[k] = d_speed;
  }

  void get_upwinded_values(double t, const mc::Vertex &v, std::vector<double> &A, std::vector<double> &Q) override {
    assert(v.get_edge_neighbors().size() == 1);
    assert(A.size() == 1);
    assert(Q.size() == 1);

    A[0] = 1.;
    Q[0] = d_speed;
  }

private:
  double d_speed;
};

class UpwindProviderNonlinearFlow : public mc::UpwindProvider {
public:
  explicit UpwindProviderNonlinearFlow(std::shared_ptr<mc::FlowUpwindEvaluator> evaluator, std::shared_ptr<mc::ExplicitNonlinearFlowSolver> solver)
      : d_evaluator(std::move(evaluator)),
        d_solver(std::move(solver)) {}

  ~UpwindProviderNonlinearFlow() override = default;

  void get_values_at_qp(double t,
                        const mc::Edge &edge,
                        size_t micro_edge,
                        const mc::QuadratureFormula &qf,
                        std::vector<double> &v_qp) override {
    assert(v_qp.size() == qf.size());

    mc::FETypeNetwork fe(qf, d_solver->get_degree());
    auto& ldof_map = d_solver->get_dof_map().get_local_dof_map(edge);
    std::vector< size_t > dof_indices(ldof_map.num_basis_functions());
    std::vector< double > dof_values(ldof_map.num_basis_functions());

    std::vector< double > values_A(qf.size());
    std::vector< double > values_Q(qf.size());

    ldof_map.dof_indices(micro_edge, d_solver->A_component, dof_indices);
    mc::extract_dof(dof_indices, d_solver->get_solution(), dof_values);
    fe.evaluate_dof_at_quadrature_points(dof_values, values_A);

    ldof_map.dof_indices(micro_edge, d_solver->Q_component, dof_indices);
    mc::extract_dof(dof_indices, d_solver->get_solution(), dof_values);
    fe.evaluate_dof_at_quadrature_points(dof_values, values_Q);

    for (size_t k = 0; k < qf.size(); k += 1)
      v_qp[k] = values_Q[k] / values_A[k];
  }

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const mc::Edge &edge, std::vector<double> &v_qp) override {
    std::vector<double> Q_up(v_qp.size());
    std::vector<double> A_up(v_qp.size());
    d_evaluator->get_fluxes_on_macro_edge(t, edge, d_solver->get_solution(), Q_up, A_up);
    for (size_t k = 0; k < v_qp.size(); k += 1)
      v_qp[k] = Q_up[k] / A_up[k];
  }

  void get_upwinded_values(double t, const mc::Vertex &v, std::vector<double> &A, std::vector<double> &Q) override {
    d_evaluator->get_fluxes_on_nfurcation(t, v, Q, A);
  }

private:
  std::shared_ptr<mc::FlowUpwindEvaluator> d_evaluator;
  std::shared_ptr<mc::ExplicitNonlinearFlowSolver> d_solver;
};

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));

  {
    const double t_end = 2.;
    const std::size_t max_iter = 1600000;
    // const std::size_t max_iter = 1;

    const double tau = 1e-4;
    const double tau_out = 1e-2;
    // const double tau_out = tau;
    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    const std::size_t num_micro_edges = 20;

    // vessel parameters
    const double vessel_length = 21.1 * 4;
    const double radius = 0.403;
    const double wall_thickness = 0.067;
    const double elastic_modulus = 400000.0;
    const double density = 1.028e-3;

    auto physical_data = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 9., radius, vessel_length);

    std::cout << mc::calculate_c0(physical_data.G0, physical_data.rho, physical_data.A0) << std::endl;

    // create_for_node the geometry of the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();
    auto v0 = graph->create_vertex();
    auto v1 = graph->create_vertex();

    auto edge_0 = graph->connect(*v0, *v1, num_micro_edges);
    edge_0->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(1., 0, 0)}});
    edge_0->add_physical_data(physical_data);

    // v0->set_to_inflow([](double t) { return 1.; });
    v0->set_to_inflow(mc::heart_beat_inflow(4.));
    v1->set_to_free_outflow();

    graph->finalize_bcs();

    // partition graph
    // TODO: app crashes if not enabled -> fix this!
    mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

    // configure solver
    auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);

    auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

    auto flow_solver = std::make_shared< mc::ExplicitNonlinearFlowSolver > (MPI_COMM_WORLD, graph, dof_map_flow, degree);
    flow_solver->use_ssp_method();

    auto upwind_evaluator = std::make_shared< mc::FlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
    auto upwind_provider = std::make_shared<UpwindProviderNonlinearFlow>(upwind_evaluator, flow_solver);
    // auto upwind_provider = std::make_shared<ConstantUpwindProvider>(vessel_length);

    mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, upwind_provider, degree);

    mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

    const auto begin_t = std::chrono::steady_clock::now();
    double t = 0;
    for (std::size_t it = 0; it < max_iter; it += 1) {

      flow_solver->solve(tau, t);
      transport_solver.solve(tau, t);

      t += tau;

      if (it % output_interval == 0) {
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "iter = " << it << ", time = " << t << std::endl;

        // save solution
        std::vector<mc::Point> points;
        std::vector<double> c_vertex_values;
        std::vector<double> Q_vertex_values;
        std::vector<double> A_vertex_values;
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->A_component, flow_solver->get_solution(), points, A_vertex_values);
        mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->Q_component, flow_solver->get_solution(), points, Q_vertex_values);

        pvd_writer.set_points(points);
        pvd_writer.add_vertex_data("c", c_vertex_values);
        pvd_writer.add_vertex_data("Q", Q_vertex_values);
        pvd_writer.add_vertex_data("A", A_vertex_values);
        pvd_writer.write(t);
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
