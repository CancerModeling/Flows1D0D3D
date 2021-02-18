////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/libmesh.h"
#include <cmath>
#include <graph_data_writer.hpp>
#include <memory>

#include "../systems/explicit_nonlinear_flow_solver.hpp"
#include "../systems/graph_storage.hpp"
#include "../systems/right_hand_side_evaluator.hpp"
#include "../systems/vessel_data_storage.hpp"
#include "../systems/vessel_formulas.hpp"
#include "../systems/errornorm.hpp"

namespace lm = libMesh;
namespace mc = macrocirculation;

class test_S {
public:
  test_S(double length, double c0, double A0)
      : d_length(length), d_c0(c0), d_A0(A0) {}

  void operator()(double t, const std::vector<lm::Point> &ps, const std::vector<double> &Q, const std::vector<double> &A, std::vector<double> &S_Q_out, std::vector<double> &S_A_out) const {
    // all vectors have to have the same shape
    assert(Q.size() == A.size());
    assert(Q.size() == S_Q_out.size());
    assert(Q.size() == S_A_out.size());
    assert(Q.size() == ps.size());

    for (std::size_t qp = 0; qp < Q.size(); qp += 1) {
      const double x = ps[qp](0);
      S_Q_out[qp] = 2 * M_PI * std::pow(d_c0, 2) * t * std::cos(M_PI * x / d_length) * std::exp(-10. * t) * std::pow(t * std::sin(M_PI * x / d_length) * std::exp(-10. * t) + 1, 2) / (std::sqrt(d_A0) * d_length);
      S_A_out[qp] = -2 * std::sin(M_PI * x / d_length) * std::exp(-10. * t) * (10. * t - 1) * (t * std::sin(M_PI * x / d_length) * std::exp(-10. * t) + 1);
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


template <std::size_t degree>
double run_scenario(std::size_t m) {
  const double t_end = 0.1;
  const std::size_t max_iter = 1000000;

  const double K_CFL = 0.25;

  const std::size_t num_edges_per_segment = 1 << (3 + m);

  const double G0 = 15.41e+1;
  const double A0 = 1.;
  const double rho = 1.028;
  const double c0 = std::sqrt(G0 / (2 * rho));
  const double length = 20.;

  const double h = length / num_edges_per_segment;

  const double tau = get_tau(h, K_CFL, c0);
  const auto output_interval = 1000000;

  std::cout << "tau " << tau << " " << c0 << std::endl;

  // we create data for the ascending aorta
  auto vessel_data = std::make_shared<mc::VesselDataStorage>();
  std::size_t ascending_aorta_id = vessel_data->add_parameter({G0, A0, rho});

  // create the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto start = graph->create_vertex(lm::Point(0, 0, 0));
  auto end = graph->create_vertex(lm::Point(length, 0, 0));
  graph->line_to(*start, *end, ascending_aorta_id, num_edges_per_segment);

  // set inflow boundary conditions
  start->set_to_inflow([](auto) { return 0.; });
  end->set_to_inflow([](auto) { return 0.; });

  // configure solver
  mc::ExplicitNonlinearFlowSolver<degree> solver(graph, vessel_data);
  solver.set_tau(tau);
  solver.use_ssp_method();
  solver.get_rhs_evaluator().set_rhs_S(test_S(length, c0, A0));

  std::vector<double> Q_vertex_values(graph->num_edges() * 2, 0);
  std::vector<double> A_vertex_values(graph->num_edges() * 2, 0);

  for (std::size_t it = 0; it < max_iter; it += 1) {
    std::cout << "iter " << it << std::endl;

    solver.solve();

    if (it % output_interval == 0) {
      // save solution
      solver.get_solution_on_vertices(Q_vertex_values, A_vertex_values);
      mc::GraphDataWriter writer;
      writer.add_vertex_data("Q", Q_vertex_values);
      writer.add_vertex_data("A", A_vertex_values);
      writer.write_vtk("line_convergence", *graph, it);
    }

    // break
    if (solver.get_time() > t_end + 1e-12)
      break;
  }

  const double error_Q = mc::errornorm<degree>(*graph, solver.get_dof_map(), 1, solver.get_solution(), [&solver,length](const std::vector< lm::Point > & p, std::vector< double >& out){
    for (std::size_t qp=0; qp<p.size(); qp+=1)
      out[qp] = get_analytic_solution_A(solver.get_time(), p[qp](0), length);
  });
  const double error_A = mc::errornorm<degree>(*graph, solver.get_dof_map(), 0, solver.get_solution(), [&solver,length](const std::vector< lm::Point > & p, std::vector< double >& out){
    for (std::size_t qp=0; qp<p.size(); qp+=1)
      out[qp] = 0;
  });

  const double error = std::sqrt(std::pow(error_Q, 2) + std::pow(error_A, 2));

  return error;
}

template < std::size_t degree >
void run_convergence_study(std::size_t max_m)
{
  std::ofstream f("convergence_error_dg" + std::to_string(degree) +  ".csv");
  for (std::size_t m=0; m<max_m; m+=1)
  {
    const double error = run_scenario<degree>(m);
    const double h = 20. / (1<<(3+m));
    f << h << ", " << error << std::endl;
  }
}

int main(int argc, char *argv[]) {
  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  lm::LibMeshInit init(argc, argv);

  const std::size_t max_m = 8;

  run_convergence_study<0>(max_m);
  run_convergence_study<1>(max_m);
  run_convergence_study<2>(max_m);
  run_convergence_study<3>(max_m);
}
