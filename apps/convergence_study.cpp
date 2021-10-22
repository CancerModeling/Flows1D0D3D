////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <macrocirculation/interpolate_to_vertices.hpp>
#include <memory>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/errornorm.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/right_hand_side_evaluator.hpp"
#include "macrocirculation/vessel_formulas.hpp"

namespace mc = macrocirculation;

class test_S {
public:
  test_S(double length, double c0, double A0)
      : d_length(length),
        d_c0(c0),
        d_A0(A0) {}

  void operator()(double t,
                  const mc::Edge &,
                  const std::vector<double> &ps,
                  const std::vector<double> &Q,
                  const std::vector<double> &A,
                  std::vector<double> &S_Q_out,
                  std::vector<double> &S_A_out) const {
    // all vectors have to have the same shape
    assert(Q.size() == A.size());
    assert(Q.size() == S_Q_out.size());
    assert(Q.size() == S_A_out.size());
    assert(Q.size() == ps.size());

    for (std::size_t qp = 0; qp < Q.size(); qp += 1) {
      const double x = ps[qp];
      S_Q_out[qp] = 2 * M_PI * std::pow(d_c0, 2) * t * std::cos(M_PI * x / d_length) * std::exp(-10. * t) *
                    std::pow(t * std::sin(M_PI * x / d_length) * std::exp(-10. * t) + 1, 2) /
                    (std::sqrt(d_A0) * d_length);
      S_A_out[qp] = -2 * std::sin(M_PI * x / d_length) * std::exp(-10. * t) * (10. * t - 1) *
                    (t * std::sin(M_PI * x / d_length) * std::exp(-10. * t) + 1);
    }
  }

private:
  double d_length;
  double d_c0;
  double d_A0;
};

double get_tau(double h, double K_CFL, double c0) {
  return h * K_CFL / c0 * 0.5;
}

double get_analytic_solution_A(double t, double x, double length) {
  return std::pow(1 + t * std::exp(-10 * t) * std::sin(M_PI * x / length), 2);
}

const std::size_t max_iter = 1000000;
const double G0 = 15.41e+1;
const double A0 = 1.;
const double rho = 1.028;
const double c0 = std::sqrt(G0 / (2 * rho));
const double length = 20.;
const double K_CFL = 0.25;
const double t_end = 0.1;

template<std::size_t degree>
double run_scenario(std::size_t num_micro_edges_per_segment, double tau, bool use_ssp) {
  // we create_for_node data for the ascending aorta
  // auto vessel_data = std::make_shared<mc::VesselDataStorage>();
  // std::size_t ascending_aorta_id = vessel_data->add_parameter({G0, A0, rho});

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto start = graph->create_vertex();
  auto end = graph->create_vertex();
  auto vessel = graph->connect(*start, *end, num_micro_edges_per_segment);
  vessel->add_embedding_data(mc::EmbeddingData{{mc::Point(0, 0, 0), mc::Point(length, 0, 0)}});
  vessel->add_physical_data(mc::PhysicalData{0., G0, A0, rho, length, 4.5e-2, 9, std::sqrt(A0 / M_PI)});

  // set inflow boundary conditions
  start->set_to_inflow_with_fixed_flow([](auto) { return 0.; });
  end->set_to_inflow_with_fixed_flow([](auto) { return 0.; });

  graph->finalize_bcs();

  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  //mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  // configure solver
  mc::ExplicitNonlinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map, degree);
  if (use_ssp)
    solver.use_ssp_method();
  else
    solver.use_explicit_euler_method();
  solver.get_rhs_evaluator().set_rhs_S(test_S(length, c0, A0));

  std::vector<double> Q_vertex_values;
  std::vector<double> A_vertex_values;
  std::vector<mc::Point> points;

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "spatial_convergence_study");

  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {
    solver.solve(tau, t);
    t += tau;

    mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 0, solver.get_solution(), points, Q_vertex_values);
    mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map, 1, solver.get_solution(), points, A_vertex_values);

    pvd_writer.set_points(points);
    pvd_writer.add_vertex_data("Q", Q_vertex_values);
    pvd_writer.add_vertex_data("A", A_vertex_values);
    pvd_writer.write(t);

    // break
    if (t >= t_end - 1e-12)
      break;
  }

  const double error_A = mc::errornorm(MPI_COMM_WORLD,
                                       *graph,
                                       solver.get_dof_map(),
                                       1,
                                       solver.get_solution(),
                                       [&solver, &t](const std::vector<double> &p, std::vector<double> &out) {
                                         for (std::size_t qp = 0; qp < p.size(); qp += 1)
                                           out[qp] = get_analytic_solution_A(t, p[qp], length);
                                       });
  const double error_Q = mc::errornorm(MPI_COMM_WORLD,
                                       *graph,
                                       solver.get_dof_map(),
                                       0,
                                       solver.get_solution(),
                                       [](const std::vector<double> &p, std::vector<double> &out) {
                                         for (std::size_t qp = 0; qp < p.size(); qp += 1)
                                           out[qp] = 0;
                                       });

  return std::sqrt(std::pow(error_Q, 2) + std::pow(error_A, 2));
}

template<std::size_t degree>
void run_temporal_convergence_study(std::size_t num_micro_edges_per_segment, std::size_t m_max) {
  std::cout << "run temporal convergence study" << std::endl;

  // we create_for_node data for the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto start = graph->create_vertex();
  auto end = graph->create_vertex();
  auto vessel = graph->connect(*start, *end, num_micro_edges_per_segment);
  vessel->add_embedding_data(mc::EmbeddingData{{mc::Point(0, 0, 0), mc::Point(length, 0, 0)}});
  vessel->add_physical_data(mc::PhysicalData{0, G0, A0, rho, length, 4.5e-2, 9, std::sqrt(A0 / M_PI)});

  // set inflow boundary conditions
  start->set_to_inflow_with_fixed_flow([](auto) { return 0.; });
  end->set_to_inflow_with_fixed_flow([](auto) { return 0.; });

  graph->finalize_bcs();

  auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  //mc::naive_mesh_partitioner( *graph, MPI_COMM_WORLD );

  const double h = 20. / num_micro_edges_per_segment;
  const double tau_max = get_tau(h, K_CFL, c0);
  const double tau_min = tau_max / (1 << (m_max + 2));
  const double t_end = std::ceil(0.1 / tau_max) * tau_max;

  // get reference solution for smallest time step.
  std::vector<double> reference_solution;
  {
    mc::ExplicitNonlinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map, degree);
    solver.get_rhs_evaluator().set_rhs_S(test_S(length, c0, A0));
    solver.use_ssp_method();

    double t = 0;
    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver.solve(tau_min, t);
      t += tau_min;
      if (t >= t_end - 1e-12)
        break;
    }
    reference_solution = solver.get_solution();
    std::cout << "finished calculating reference solution" << std::endl;
  }

  std::vector<double> diff(reference_solution.size(), 0);

  const auto zero_fct = [](const std::vector<double> &p, std::vector<double> &out) {
    for (std::size_t qp = 0; qp < p.size(); qp += 1)
      out[qp] = 0;
  };

  std::ofstream f("temporal_convergence_error_dg" + std::to_string(degree) + ".csv");

  for (std::size_t m = 0; m < m_max; m += 1) {
    const auto begin_t = std::chrono::steady_clock::now();

    const double tau = tau_max / (1 << m);

    mc::ExplicitNonlinearFlowSolver solver(MPI_COMM_WORLD, graph, dof_map, degree);
    solver.get_rhs_evaluator().set_rhs_S(test_S(length, c0, A0));
    solver.use_ssp_method();

    double t = 0;
    for (std::size_t it = 0; it < max_iter; it += 1) {
      solver.solve(tau, t);
      t += tau;

      if (t >= t_end - 1e-12)
        break;
    }

    gmm::add(reference_solution, gmm::scaled(solver.get_solution(), -1), diff);

    const double error_Q = mc::errornorm(MPI_COMM_WORLD, *graph, solver.get_dof_map(), 0, diff, zero_fct);
    const double error_A = mc::errornorm(MPI_COMM_WORLD, *graph, solver.get_dof_map(), 1, diff, zero_fct);
    const double error = std::sqrt(std::pow(error_Q, 2) + std::pow(error_A, 2));

    if (mc::mpi::rank(MPI_COMM_WORLD) == 0) {
      f << tau << ", " << error << std::endl;

      const auto end_t = std::chrono::steady_clock::now();
      const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
      std::cout << "finished tau = " << tau << " (error = " << error << ", time = " << elapsed_ms * 1e-6 << " s)" << std::endl;
    }
  }
}

template<std::size_t degree>
void run_spatial_convergence_study(std::size_t max_m) {
  std::cout << "spatial convergence study for degree " << std::to_string(degree) << std::endl;
  std::ofstream f("spatial_convergence_error_dg" + std::to_string(degree) + ".csv");
  const double tau = get_tau(20. / (1 << (3 + max_m)), K_CFL, c0);
  double previous_error = NAN;
  for (std::size_t m = 0; m < max_m; m += 1) {
    const auto begin = std::chrono::steady_clock::now();

    const std::size_t num_edges_per_segment = 1 << (3 + m);
    const double error = run_scenario<degree>(num_edges_per_segment, tau, true);
    const double h = 20. / num_edges_per_segment;
    const double rate = error / previous_error;
    f << h << ", " << error << std::endl;

    const auto end = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    std::cout << "finished h = " << h << " (error = " << error << ", time = " << elapsed_ms * 1e-6 << " s, rate = " << rate << ")" << std::endl;

    previous_error = error;
  }
}

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  MPI_Init(&argc, &argv);

  std::cout << "comm world size = " << mc::mpi::size(MPI_COMM_WORLD) << std::endl;

  //const std::size_t max_m = 8;
  const std::size_t max_m = 8;
  run_spatial_convergence_study<0>(max_m);
  run_spatial_convergence_study<1>(max_m);
  run_spatial_convergence_study<2>(max_m);
  run_spatial_convergence_study<3>(max_m);

  run_temporal_convergence_study<0>(12, 8);
  run_temporal_convergence_study<1>(12, 8);
  run_temporal_convergence_study<2>(12, 8);
  run_temporal_convergence_study<3>(12, 8);

  MPI_Finalize();
}
