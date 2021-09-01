////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <memory>
#include <utility>

#include "petsc.h"

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/explicit_transport_solver.hpp"
#include "macrocirculation/fe_type.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/right_hand_side_evaluator.hpp"
#include "macrocirculation/vessel_formulas.hpp"


namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

namespace linearized {

/*! @brief Calculates the upwinding at an inner boundary.
 *
 * @param alpha The \f$ \sqrt{\frac{C}{L}} \f$ factor for the linearized characteristics.
 * @param p_l   The pressure at the boundary of the left cell.
 * @param q_l   The flux at the boundary of the left cell.
 * @param p_r   The pressure at the boundary of the right cell.
 * @param q_r   The flux at the boundary of the right cell.
 * @param p_up  The upwinded pressure.
 * @param q_up  The upwinded flux.
 */
inline void inner_boundary(double alpha, double p_l, double q_l, double p_r, double q_r, double &p_up, double &q_up) {
  p_up = 0.5 * (p_l + p_r + alpha * (q_l - q_r));
  q_up = alpha * (-p_up + p_l) + q_l;
}

/*! @brief Calculates the upwinding at an nfurcation.
 *
 * @param p         The boundary pressure values of each vessel at a common vertex.
 * @param q         The boundary flux values of each vessel at a common vertex.
 * @param params    The physical edge parameters near the nfurcation.
 * @param sigma     The vessel normals (+1 if the vessel points to the vertex, -1 if not).
 * @param p_up      The upwinded pressure values.
 * @param q_up      The upwinded flux values.
 */
inline void nfurcation_boundary(const std::vector<double> &p,
                                const std::vector<double> &q,
                                const std::vector<mc::VesselParameters> &params,
                                const std::vector<double> &sigma,
                                std::vector<double> &p_up,
                                std::vector<double> &q_up) {
  // all vectors need to have the same size:
  assert(p.size() == q.size());
  assert(p_up.size() == p.size());
  assert(q_up.size() == p.size());
  assert(param.size() == p.size());
  assert(in.size() == p.size());

  const size_t N = p.size();

  // dof-ordering: (p_1, q_1, p_2, q_2, ... p_N, q_N)
  Eigen::MatrixXd mat(p.size() + q.size(), p.size() + q.size());
  Eigen::VectorXd rhs(p.size() + q.size());

  mat.setZero();

  // vector containing the \f$ \sqrt{ \frac{C}{L} } \f$ factors:
  std::vector<double> alpha;
  for (auto &param : params) {
    alpha.push_back(std::sqrt(mc::linear::get_C(param) / mc::linear::get_L(param)));
  }

  // constrain the upwinded fluxes to zero:
  for (int k = 0; k < N; k += 1) {
    mat(0, 2 * k + 1) = sigma[k];
  }
  rhs(0) = 0;

  // the characteristics should be equal
  for (int k = 0; k < N; k += 1) {
    mat(1 + k, 2 * k) = 0.5 * alpha[k] * sigma[k];
    mat(1 + k, 2 * k + 1) = 0.5;
    rhs(1 + k) = 0.5 * alpha[k] * sigma[k] * p[k] + 0.5 * q[k];
  }

  // the pressures should be equal
  for (int k = 1; k < N; k += 1) {
    mat(N + k, 2 * k) = 1.;
    mat(N + k, 2 * (k - 1)) = -1.;
    rhs(N + k) = 0.;
  }

  // std::cout << mat << std::endl;
  // std::cout << rhs << std::endl;

  Eigen::VectorXd result = mat.fullPivLu().solve(rhs);

  for (int k = 0; k < N; k += 1) {
    p_up[k] = result[2 * k];
    q_up[k] = result[2 * k + 1];
  }
}

} // namespace linearized

template<typename VectorType>
inline void compute_inner_boundary_values_on_macro_edge(const mc::Edge &edge,
                                                        const mc::DofMap &dof_map,
                                                        const VectorType &u,
                                                        size_t component_index,
                                                        std::vector<double> &boundary_values_l,
                                                        std::vector<double> &boundary_values_r) {
  auto &local_dof_map = dof_map.get_local_dof_map(edge);

  assert(local_dof_map.num_micro_edges() == boundary_values.size());

  const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

  const auto &param = edge.get_physical_data();
  const double h = param.length / local_dof_map.num_micro_edges();

  std::vector<std::size_t> dof_indices(num_basis_functions, 0);
  std::vector<double> dof_values(num_basis_functions, 0);

  mc::FETypeNetwork fe(mc::create_trapezoidal_rule(), num_basis_functions - 1);
  fe.reinit(h);

  for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_vertices(); micro_edge_id += 1) {
    local_dof_map.dof_indices(micro_edge_id, component_index, dof_indices);
    macrocirculation::extract_dof(dof_indices, u, dof_values);
    auto boundary_values = fe.evaluate_dof_at_boundary_points(dof_values);
    boundary_values_l[micro_edge_id] = boundary_values.left;
    boundary_values_r[micro_edge_id] = boundary_values.right;
  }
}

class LinearizedFlowUpwindEvaluator {
public:
  LinearizedFlowUpwindEvaluator(MPI_Comm comm, std::shared_ptr<mc::GraphStorage> graph, std::shared_ptr<mc::DofMap> dof_map)
      : d_comm(comm),
        d_graph(std::move(graph)),
        d_dof_map(std::move(dof_map)),
        d_q_boundary_evaluator(comm, d_graph, d_dof_map, 1),
        d_p_boundary_evaluator(comm, d_graph, d_dof_map, 0),
        d_q_macro_edge_flux_l(d_graph->num_edges()),
        d_q_macro_edge_flux_r(d_graph->num_edges()),
        d_p_macro_edge_flux_l(d_graph->num_edges()),
        d_p_macro_edge_flux_r(d_graph->num_edges()),
        d_current_t(NAN) {}

  void init(double t, const std::vector<double> &u_prev) {
    init_generic(t, u_prev);
  }

  void init(double t, const mc::PetscVec &u_prev) {
    init_generic(t, u_prev);
  }

  void get_fluxes_on_macro_edge(double t, const mc::Edge &edge, const std::vector<double> &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
    get_fluxes_on_macro_edge_generic(t, edge, u_prev, p_up, q_up);
  }

  void get_fluxes_on_macro_edge(double t, const mc::Edge &edge, const mc::PetscVec &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
    get_fluxes_on_macro_edge_generic(t, edge, u_prev, p_up, q_up);
  }

  void get_fluxes_on_nfurcation(double t, const mc::Vertex &v, std::vector<double> &p_up, std::vector<double> &q_up) const {
    // evaluator was initialized with the correct time step
    if (d_current_t != t)
      throw std::runtime_error("LinearizedFlowUpwindEvaluator was not initialized for the given time step");

    p_up.resize(v.get_edge_neighbors().size());
    q_up.resize(v.get_edge_neighbors().size());

    for (size_t neighbor_edge_idx = 0; neighbor_edge_idx < v.get_edge_neighbors().size(); neighbor_edge_idx += 1) {
      const auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[neighbor_edge_idx]);

      if (edge.is_pointing_to(v.get_id())) {
        p_up[neighbor_edge_idx] = d_p_macro_edge_flux_r[edge.get_id()];
        q_up[neighbor_edge_idx] = d_q_macro_edge_flux_r[edge.get_id()];
      } else {
        p_up[neighbor_edge_idx] = d_p_macro_edge_flux_l[edge.get_id()];
        q_up[neighbor_edge_idx] = d_q_macro_edge_flux_l[edge.get_id()];
      }
    }
  }

private:
  template<typename VectorType>
  void init_generic(double t, const VectorType &u_prev) {
    d_current_t = t;

    d_p_boundary_evaluator.init(u_prev);
    d_q_boundary_evaluator.init(u_prev);

    calculate_nfurcation_fluxes(u_prev);
    calculate_inout_fluxes(t, u_prev);
  }

  template<typename VectorType>
  void get_fluxes_on_macro_edge_generic(double t, const mc::Edge &edge, const VectorType &u_prev, std::vector<double> &p_up, std::vector<double> &q_up) const {
    // evaluator was initialized with the correct time step
    if (d_current_t != t)
      throw std::runtime_error("LinearizedFlowUpwindEvaluator was not initialized for the given time step");

    assert(p_up_macro_edge.size() == edge.num_micro_vertices());
    assert(q_up_macro_edge.size() == edge.num_micro_vertices());

    std::vector<double> p_l(edge.num_micro_edges(), 0);
    std::vector<double> p_r(edge.num_micro_edges(), 0);
    std::vector<double> q_l(edge.num_micro_edges(), 0);
    std::vector<double> q_r(edge.num_micro_edges(), 0);

    compute_inner_boundary_values_on_macro_edge(edge, *d_dof_map, u_prev, 0, p_l, p_r);
    compute_inner_boundary_values_on_macro_edge(edge, *d_dof_map, u_prev, 1, q_l, q_r);

    assert(edge.has_physical_data());
    auto param = edge.get_physical_data();

    // TODO: move this calculation into its own function
    const auto alpha = std::sqrt(mc::linear::get_C(param) / mc::linear::get_L(param));

    for (std::size_t micro_vertex_id = 1; micro_vertex_id < edge.num_micro_vertices() - 1; micro_vertex_id += 1) {
      linearized::inner_boundary(alpha, p_l[micro_vertex_id], q_l[micro_vertex_id], p_r[micro_vertex_id - 1], q_r[micro_vertex_id - 1], p_up[micro_vertex_id], q_up[micro_vertex_id]);
    }

    // update left fluxes
    p_up[0] = d_p_macro_edge_flux_l[edge.get_id()];
    q_up[0] = d_q_macro_edge_flux_l[edge.get_id()];

    // update right fluxes
    p_up[edge.num_micro_vertices() - 1] = d_p_macro_edge_flux_r[edge.get_id()];
    q_up[edge.num_micro_vertices() - 1] = d_q_macro_edge_flux_r[edge.get_id()];
  }

  template<typename VectorType>
  void calculate_nfurcation_fluxes(const VectorType &u_prev) {
    for (const auto &v_id : d_graph->get_active_vertex_ids(mc::mpi::rank(d_comm))) {
      const auto vertex = d_graph->get_vertex(v_id);

      // we only handle bifurcations
      if (!vertex->is_bifurcation())
        continue;

      const size_t num_vessels = vertex->get_edge_neighbors().size();

      // get edges
      std::vector<std::shared_ptr<mc::Edge>> e;
      for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
        e.push_back(d_graph->get_edge(vertex->get_edge_neighbors()[vessel_idx]));

      // check orientation
      std::vector<double> sigma;
      for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
        sigma.push_back(e[vessel_idx]->is_pointing_to(vertex->get_id()) ? +1 : -1);

      // get data
      std::vector<mc::VesselParameters> params;
      for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
        const auto &data_e = e[vessel_idx]->get_physical_data();
        params.emplace_back(data_e.G0, data_e.A0, data_e.rho);
      }

      std::vector<double> p_e(num_vessels);
      std::vector<double> q_e(num_vessels);
      d_p_boundary_evaluator(*vertex, p_e);
      d_q_boundary_evaluator(*vertex, q_e);

      std::vector<double> p_up(num_vessels, 0);
      std::vector<double> q_up(num_vessels, 0);

      linearized::nfurcation_boundary(p_e, q_e, params, sigma, p_up, q_up);

      for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
        if (e[vessel_idx]->is_pointing_to(v_id)) {
          d_p_macro_edge_flux_r[e[vessel_idx]->get_id()] = p_up[vessel_idx];
          d_q_macro_edge_flux_r[e[vessel_idx]->get_id()] = q_up[vessel_idx];
        } else {
          d_p_macro_edge_flux_l[e[vessel_idx]->get_id()] = p_up[vessel_idx];
          d_q_macro_edge_flux_l[e[vessel_idx]->get_id()] = q_up[vessel_idx];
        }
      }
    }
  }

  template<typename VectorType>
  void calculate_inout_fluxes(double t, const VectorType &u_prev) {
    for (const auto &v_id : d_graph->get_active_vertex_ids(mc::mpi::rank(d_comm))) {
      const auto vertex = d_graph->get_vertex(v_id);

      // we are only interested in leaves and ignore the rest:
      if (!vertex->is_leaf())
        continue;

      const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);
      const auto &param = edge->get_physical_data();
      const double sigma = edge->is_pointing_to(vertex->get_id()) ? +1. : -1.;

      std::vector<double> p_value;
      std::vector<double> q_value;
      d_p_boundary_evaluator(*vertex, p_value);
      d_q_boundary_evaluator(*vertex, q_value);

      // TODO: move this calculation into its own function
      const auto alpha = std::sqrt(mc::linear::get_C(param) / mc::linear::get_L(param));

      double p_up = NAN;
      double q_up = NAN;

      if (vertex->is_inflow()) {
        q_up = -sigma * vertex->get_inflow_value(t);
        p_up = p_value[0] + 1. / alpha * (q_value[0] - q_up);
      } else if (vertex->is_free_outflow()) {
        p_up = 0.5 * (p_value[0] + sigma / alpha * q_value[0]);
        q_up = sigma * alpha * p_up;
      } else {
        throw std::runtime_error("boundary type not implemented yet in LinearizedFlowUpwindEvaluator.");
      }

      // write back into vectors:
      if (edge->is_pointing_to(vertex->get_id())) {
        d_q_macro_edge_flux_r[edge->get_id()] = q_up;
        d_p_macro_edge_flux_r[edge->get_id()] = p_up;
      } else {
        d_q_macro_edge_flux_l[edge->get_id()] = q_up;
        d_p_macro_edge_flux_l[edge->get_id()] = p_up;
      }
    }
  }

private:
  MPI_Comm d_comm;

  std::shared_ptr<mc::GraphStorage> d_graph;

  std::shared_ptr<mc::DofMap> d_dof_map;

  mc::EdgeBoundaryEvaluator d_p_boundary_evaluator;

  mc::EdgeBoundaryEvaluator d_q_boundary_evaluator;

  std::vector<double> d_p_macro_edge_flux_l;
  std::vector<double> d_p_macro_edge_flux_r;
  std::vector<double> d_q_macro_edge_flux_l;
  std::vector<double> d_q_macro_edge_flux_r;

  double d_current_t;
};

class UpwindProviderLinearizedFlow : public mc::UpwindProvider {
public:
  explicit UpwindProviderLinearizedFlow(std::shared_ptr<mc::GraphStorage> graph, std::shared_ptr<LinearizedFlowUpwindEvaluator> evaluator, std::shared_ptr<mc::ImplicitLinearFlowSolver> solver)
      : d_graph(std::move(graph)),
        d_evaluator(std::move(evaluator)),
        d_solver(std::move(solver)) {}

  ~UpwindProviderLinearizedFlow() override = default;

  void init(double t, const std::vector<double> &u) override {
    d_evaluator->init(t, u);
  }

  void init(double t, const mc::PetscVec &u) override {
    d_evaluator->init(t, u);
  }

  void get_values_at_qp(double t,
                        const mc::Edge &edge,
                        size_t micro_edge,
                        const mc::QuadratureFormula &qf,
                        std::vector<double> &v_qp) const override {
    assert(v_qp.size() == qf.size());

    mc::FETypeNetwork fe(qf, d_solver->get_degree());
    auto &ldof_map = d_solver->get_dof_map().get_local_dof_map(edge);
    std::vector<size_t> dof_indices(ldof_map.num_basis_functions());
    std::vector<double> dof_values(ldof_map.num_basis_functions());

    std::vector<double> values_q(qf.size());

    ldof_map.dof_indices(micro_edge, d_solver->q_component, dof_indices);
    mc::extract_dof(dof_indices, d_solver->get_solution(), dof_values);
    fe.evaluate_dof_at_quadrature_points(dof_values, values_q);

    auto &param = edge.get_physical_data();

    for (size_t k = 0; k < qf.size(); k += 1)
      v_qp[k] = values_q[k] / param.A0;

    // std::cout << v_qp << std::endl;
  }

  /*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
  void get_upwinded_values(double t, const mc::Edge &edge, std::vector<double> &v_qp) const override {
    assert(v_qp.size() == edge.num_micro_vertices());
    std::vector<double> p_up(edge.num_micro_vertices());
    std::vector<double> q_up(edge.num_micro_vertices());
    d_evaluator->get_fluxes_on_macro_edge(t, edge, d_solver->get_solution(), p_up, q_up);
    assert(edge.has_physical_data());
    auto A0 = edge.get_physical_data().A0;
    for (size_t k = 0; k < v_qp.size(); k += 1)
      v_qp[k] = q_up[k] / A0;

    // std::cout << v_qp << std::endl;
  }

  void get_upwinded_values(double t, const mc::Vertex &v, std::vector<double> &A, std::vector<double> &Q) const override {
    std::vector<double> p_up(v.get_edge_neighbors().size());
    d_evaluator->get_fluxes_on_nfurcation(t, v, p_up, Q);

    std::vector<double> A0;
    for (size_t k = 0; k < v.get_edge_neighbors().size(); k += 1) {
      auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[k]);
      assert(edge.has_physical_data());
      A[k] = edge.get_physical_data().A0;
    }

    // std::cout << "v = " << v.get_id() << " " << A << " " << Q << std::endl;
  }

private:
  std::shared_ptr<mc::GraphStorage> d_graph;

  std::shared_ptr<LinearizedFlowUpwindEvaluator> d_evaluator;

  std::shared_ptr<mc::ImplicitLinearFlowSolver> d_solver;
};

void implicit_transport_with_implicit_flow(double tau, double tau_out, double t_end, std::shared_ptr<mc::GraphStorage> graph) {
  const std::size_t max_iter = 1600000;

  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  // dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, true);

  auto flow_solver = std::make_shared<mc::ImplicitLinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver->setup(tau);

  auto upwind_evaluator = std::make_shared<LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
  auto variable_upwind_provider = std::make_shared<UpwindProviderLinearizedFlow>(graph, upwind_evaluator, flow_solver);

  mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    flow_solver->solve(tau, t + tau);
    variable_upwind_provider->init(t + tau, flow_solver->get_solution());
    transport_solver.solve(tau, t + tau);

    t += tau;

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", time = " << t << std::endl;

      // save solution
      std::vector<mc::Point> points;
      std::vector<double> c_vertex_values;
      std::vector<double> p_vertex_values;
      std::vector<double> q_vertex_values;
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_transport, 0, transport_solver.get_solution(), points, c_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->p_component, flow_solver->get_solution(), points, p_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->q_component, flow_solver->get_solution(), points, q_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.add_vertex_data("q", q_vertex_values);
      pvd_writer.add_vertex_data("p", p_vertex_values);
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

void implicit_transport_with_explicit_flow(double tau, double tau_out, double t_end, std::shared_ptr<mc::GraphStorage> graph) {
  const std::size_t max_iter = 1600000;

  const auto output_interval = static_cast<std::size_t>(tau_out / tau);

  // configure solver
  auto dof_map_transport = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_transport->create(MPI_COMM_WORLD, *graph, 1, degree, true);

  auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

  auto flow_solver = std::make_shared<mc::ExplicitNonlinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
  flow_solver->use_ssp_method();

  auto upwind_evaluator = std::make_shared<mc::FlowAQUpwindEvaluator>(MPI_COMM_WORLD, graph, dof_map_flow);
  auto variable_upwind_provider = std::make_shared<mc::UpwindProviderNonlinearFlow>(upwind_evaluator, flow_solver);

  mc::ImplicitTransportSolver transport_solver(MPI_COMM_WORLD, graph, dof_map_transport, variable_upwind_provider, degree);

  mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, "output", "transport_solution");

  const auto begin_t = std::chrono::steady_clock::now();
  double t = 0;
  for (std::size_t it = 0; it < max_iter; it += 1) {

    flow_solver->solve(tau, t + tau);
    variable_upwind_provider->init(t + tau, flow_solver->get_solution());
    transport_solver.solve(tau, t + tau);

    t += tau;

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", time = " << t << std::endl;

      // save solution
      std::vector<mc::Point> points;
      std::vector<double> c_vertex_values;
      std::vector<double> A_vertex_values;
      std::vector<double> Q_vertex_values;
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

int main(int argc, char *argv[]) {
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves implicit transport problem"));

  cxxopts::Options options(argv[0], "Implicit transport solver.");
  options.add_options()                                                                                                              //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                    //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("6"))                                      //
    ("no-upper-vessel", "Disables the upper vessel at the bifurcation", cxxopts::value<bool>()->default_value("false"))              //
    ("no-lower-vessel", "Disables the lower vessel at the bifurcation", cxxopts::value<bool>()->default_value("false"))              //
    ("explicit-flow", "Enables an explicit flow solver instead of the implicit one", cxxopts::value<bool>()->default_value("false")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const double tau = args["tau"].as<double>();
  const double tau_out = args["tau-out"].as<double>();
  const double t_end = args["t-end"].as<double>();

  const bool no_lower_vessel = args["no-lower-vessel"].as<bool>();
  const bool no_upper_vessel = args["no-upper-vessel"].as<bool>();

  const std::size_t num_micro_edges = 20;

  // vessel parameters
  //const double vessel_length = 20.5;
  const double vessel_length = 10.5;
  const double radius = 0.403;
  const double wall_thickness = 0.067;
  const double elastic_modulus = 400000.0;
  const double density = 1.028e-3;

  auto physical_data_short = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);
  auto physical_data_long = mc::PhysicalData::set_from_data(elastic_modulus, wall_thickness, density, 2., radius, vessel_length / 2.);

  // create_for_node the geometry of the ascending aorta
  auto graph = std::make_shared<mc::GraphStorage>();
  auto v0 = graph->create_vertex();
  auto v1 = graph->create_vertex();
  auto v2 = graph->create_vertex();

  auto edge_0 = graph->connect(*v0, *v1, num_micro_edges);
  edge_0->add_embedding_data({{mc::Point(0, 0, 0), mc::Point(0.5, 0, 0)}});
  edge_0->add_physical_data(physical_data_short);

  auto edge_1 = graph->connect(*v1, *v2, num_micro_edges);
  edge_1->add_embedding_data({{mc::Point(0.5, 0, 0), mc::Point(1., 0, 0)}});
  edge_1->add_physical_data(physical_data_long);

  if (!no_upper_vessel) {
    auto v3 = graph->create_vertex();
    auto edge_2 = graph->connect(*v3, *v1, num_micro_edges);
    edge_2->add_embedding_data({{mc::Point(0.5, 0.5, 0), mc::Point(0.5, 0, 0)}});
    edge_2->add_physical_data(physical_data_long);
    v3->set_to_free_outflow();
  }

  if (!no_lower_vessel) {
    auto v4 = graph->create_vertex();
    auto edge_3 = graph->connect(*v4, *v1, num_micro_edges);
    edge_3->add_embedding_data({{mc::Point(0.5, -0.5, 0), mc::Point(0.5, 0.0, 0)}});
    edge_3->add_physical_data(physical_data_long);
    v4->set_to_free_outflow();
  }

  v0->set_to_inflow([](double t) { return mc::heart_beat_inflow(4., 1., 0.7)(t); });
  v2->set_to_free_outflow();

  graph->finalize_bcs();

  // partition graph
  mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);

  const bool explicit_flow = args["explicit-flow"].as<bool>();

  if (explicit_flow) {
    implicit_transport_with_explicit_flow(tau, tau_out, t_end, graph);
  } else {
    implicit_transport_with_implicit_flow(tau, tau_out, t_end, graph);
  }

  CHKERRQ(PetscFinalize());
}
