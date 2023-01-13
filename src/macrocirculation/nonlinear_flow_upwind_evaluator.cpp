#include "nonlinear_flow_upwind_evaluator.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "vessel_formulas.hpp"

#include <cmath>
#include <utility>

namespace macrocirculation {

NonlinearFlowUpwindEvaluator::NonlinearFlowUpwindEvaluator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_A_boundary_evaluator(comm, d_graph, d_dof_map, 1),
      d_Q_boundary_evaluator(comm, d_graph, d_dof_map, 0),
      d_Q_macro_edge_flux_l(d_graph->num_edges()),
      d_Q_macro_edge_flux_r(d_graph->num_edges()),
      d_A_macro_edge_flux_l(d_graph->num_edges()),
      d_A_macro_edge_flux_r(d_graph->num_edges()),
      d_current_t(NAN) {}

void NonlinearFlowUpwindEvaluator::init(double t, const std::vector<double> &u_prev) {
  d_current_t = t;

  d_A_boundary_evaluator.init(u_prev);
  d_Q_boundary_evaluator.init(u_prev);

  calculate_nfurcation_fluxes(u_prev);
  calculate_inout_fluxes(t, u_prev);
}

void NonlinearFlowUpwindEvaluator::get_fluxes_on_macro_edge(double t, const Edge &edge, const std::vector<double> &u_prev, std::vector<double> &Q_up_macro_edge, std::vector<double> &A_up_macro_edge) const {
  // evaluator was initialized with the correct time step
  if (d_current_t != t)
    throw std::runtime_error("FlowUpwindEvaluator was not initialized for the given time step");

  const auto local_dof_map = d_dof_map->get_local_dof_map(edge);

  assert(Q_up_macro_edge.size() == local_dof_map.num_micro_vertices());
  assert(A_up_macro_edge.size() == local_dof_map.num_micro_vertices());

  const auto &param = edge.get_physical_data();
  const double h = param.length / local_dof_map.num_micro_edges();

  const std::size_t num_basis_functions = local_dof_map.num_basis_functions();

  // TODO: Make this more efficient by not recalculating the gradients.
  // finite-element for left and right edge
  FETypeNetwork fe(create_trapezoidal_rule(), num_basis_functions - 1);

  // dof indices for left and right edge
  std::vector<std::size_t> Q_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_l(num_basis_functions, 0);
  std::vector<std::size_t> Q_dof_indices_r(num_basis_functions, 0);
  std::vector<std::size_t> A_dof_indices_r(num_basis_functions, 0);

  // local views of our previous solution
  std::vector<double> Q_prev_loc_l(num_basis_functions, 0);
  std::vector<double> A_prev_loc_l(num_basis_functions, 0);
  std::vector<double> Q_prev_loc_r(num_basis_functions, 0);
  std::vector<double> A_prev_loc_r(num_basis_functions, 0);

  // previous solution evaluated at quadrature points
  // we have 2 quadrature points for the trapezoidal rule
  std::vector<double> Q_prev_qp_l(2, 0);
  std::vector<double> A_prev_qp_l(2, 0);
  std::vector<double> Q_prev_qp_r(2, 0);
  std::vector<double> A_prev_qp_r(2, 0);

  fe.reinit(h);

  // TODO: the boundary values of each cell are evaluated twice
  for (std::size_t micro_vertex_id = 1; micro_vertex_id < local_dof_map.num_micro_vertices() - 1; micro_vertex_id += 1) {
    const std::size_t local_micro_edge_id_l = micro_vertex_id - 1;
    const std::size_t local_micro_edge_id_r = micro_vertex_id;

    local_dof_map.dof_indices(local_micro_edge_id_l, 0, Q_dof_indices_l);
    local_dof_map.dof_indices(local_micro_edge_id_r, 0, Q_dof_indices_r);
    local_dof_map.dof_indices(local_micro_edge_id_l, 1, A_dof_indices_l);
    local_dof_map.dof_indices(local_micro_edge_id_r, 1, A_dof_indices_r);

    extract_dof(Q_dof_indices_l, u_prev, Q_prev_loc_l);
    extract_dof(A_dof_indices_l, u_prev, A_prev_loc_l);
    extract_dof(Q_dof_indices_r, u_prev, Q_prev_loc_r);
    extract_dof(A_dof_indices_r, u_prev, A_prev_loc_r);

    fe.evaluate_dof_at_quadrature_points(Q_prev_loc_l, Q_prev_qp_l);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc_l, A_prev_qp_l);
    fe.evaluate_dof_at_quadrature_points(Q_prev_loc_r, Q_prev_qp_r);
    fe.evaluate_dof_at_quadrature_points(A_prev_loc_r, A_prev_qp_r);

    const double Q_l = Q_prev_qp_l[1];
    const double A_l = A_prev_qp_l[1];
    const double Q_r = Q_prev_qp_r[0];
    const double A_r = A_prev_qp_r[0];

    const double W2_l = nonlinear::get_w2_from_QA(Q_l, A_l, param);
    const double W1_r = nonlinear::get_w1_from_QA(Q_r, A_r, param);

    double Q_up = 0, A_up = 0;

    solve_W12(Q_up, A_up, W1_r, W2_l, param.G0, param.rho, param.A0);

    Q_up_macro_edge[micro_vertex_id] = Q_up;
    A_up_macro_edge[micro_vertex_id] = A_up;
  }

  // update left fluxes
  Q_up_macro_edge[0] = d_Q_macro_edge_flux_l[edge.get_id()];
  A_up_macro_edge[0] = d_A_macro_edge_flux_l[edge.get_id()];

  // update right fluxes
  Q_up_macro_edge[local_dof_map.num_micro_vertices() - 1] = d_Q_macro_edge_flux_r[edge.get_id()];
  A_up_macro_edge[local_dof_map.num_micro_vertices() - 1] = d_A_macro_edge_flux_r[edge.get_id()];
}

void NonlinearFlowUpwindEvaluator::get_fluxes_on_nfurcation(double t, const Vertex &v, std::vector<double> &Q_up, std::vector<double> &A_up) const {
  // evaluator was initialized with the correct time step
  if (d_current_t != t)
    throw std::runtime_error("FlowUpwindEvaluator was not initialized for the given time step");

  Q_up.resize(v.get_edge_neighbors().size());
  A_up.resize(v.get_edge_neighbors().size());

  for (size_t neighbor_edge_idx = 0; neighbor_edge_idx < v.get_edge_neighbors().size(); neighbor_edge_idx += 1) {
    const auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[neighbor_edge_idx]);

    if (edge.is_pointing_to(v.get_id())) {
      Q_up[neighbor_edge_idx] = d_Q_macro_edge_flux_r[edge.get_id()];
      A_up[neighbor_edge_idx] = d_A_macro_edge_flux_r[edge.get_id()];
    } else {
      Q_up[neighbor_edge_idx] = d_Q_macro_edge_flux_l[edge.get_id()];
      A_up[neighbor_edge_idx] = d_A_macro_edge_flux_l[edge.get_id()];
    }
  }
}

void NonlinearFlowUpwindEvaluator::calculate_nfurcation_fluxes(const std::vector<double> &/*u_prev*/) {
  for (const auto &v_id : d_graph->get_active_and_connected_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // we only handle bifurcations
    if (!vertex->is_bifurcation())
      continue;

    const size_t num_vessels = vertex->get_edge_neighbors().size();

    // get edges
    std::vector<std::shared_ptr<Edge>> e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      e.push_back(d_graph->get_edge(vertex->get_edge_neighbors()[vessel_idx]));

    // check orientation
    std::vector<bool> e_in;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1)
      e_in.push_back(e[vessel_idx]->is_pointing_to(vertex->get_id()));

    // get data
    std::vector<VesselParameters> p_e;
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      const auto &data_e = e[vessel_idx]->get_physical_data();
      p_e.emplace_back(data_e.G0, data_e.A0, data_e.rho);
    }

    // evaluate on edges
    std::vector<double> Q_e(num_vessels);
    std::vector<double> A_e(num_vessels);
    d_Q_boundary_evaluator(*vertex, Q_e);
    d_A_boundary_evaluator(*vertex, A_e);

    // get upwinded values at bifurcation
    std::vector<double> Q_up(num_vessels, 0);
    std::vector<double> A_up(num_vessels, 0);
    solve_at_nfurcation(Q_e, A_e, p_e, e_in, Q_up, A_up);

    // save upwinded values into upwind vector
    for (size_t vessel_idx = 0; vessel_idx < num_vessels; vessel_idx += 1) {
      if (e_in[vessel_idx]) {
        d_Q_macro_edge_flux_r[e[vessel_idx]->get_id()] = Q_up[vessel_idx];
        d_A_macro_edge_flux_r[e[vessel_idx]->get_id()] = A_up[vessel_idx];
      } else {
        d_Q_macro_edge_flux_l[e[vessel_idx]->get_id()] = Q_up[vessel_idx];
        d_A_macro_edge_flux_l[e[vessel_idx]->get_id()] = A_up[vessel_idx];
      }
    }
  }
}

void calculate_windkessel_upwind_values(const PhysicalData &param, bool is_pointing_to, double R1, double Q_DG, double A_DG, double p_c, double &Q_out, double &A_out) {
  const auto W = is_pointing_to ? nonlinear::get_w2_from_QA(Q_DG, A_DG, param)
                                : nonlinear::get_w1_from_QA(Q_DG, A_DG, param);

  const double c0 = param.get_c0();

  const auto f = [&](auto A_out) {
    auto p = param.G0 * (std::sqrt(A_out / param.A0) - 1);
    return W - (p - p_c) / (A_out * R1) - 4 * c0 * std::pow(A_out / param.A0, 0.25);
  };

  const auto df = [&](auto A_out) {
    auto p = param.G0 * (std::sqrt(A_out / param.A0) - 1);
    auto dp = param.G0 * 0.5 / std::sqrt(A_out * param.A0);
    return -dp / (A_out * R1) + (p - p_c) / (A_out * A_out * R1) - c0 * std::pow(A_out, -0.75) / std::pow(param.A0, 0.25);
  };

  const double TOL = 1.0e-10;
  const double omega = 0.5 / 2.;
  double error = INFINITY;
  int num_iter = 0;

  A_out = A_DG;

  while (num_iter < 250 && error > TOL) {
    // std::cout << "f_prev" << f(A_out) << std::endl;
    const double f_value = f(A_out);
    const double df_value = df(A_out);
    const double dx = -f_value / df_value;

    A_out = A_out + omega * dx;

    error = std::abs(f(A_out));
    num_iter += 1;

    // std::cout << "f_now " << f(A_out) << std::endl;
    if (num_iter == 250)
      std::cerr << "warning: Newton did not converge" << std::endl;
  }

  const double sgn = is_pointing_to ? +1 : -1;

  Q_out = sgn * (nonlinear::get_p_from_A(A_out, param.G0, param.A0) - p_c) / R1;
}

void NonlinearFlowUpwindEvaluator::calculate_inout_fluxes(double t, const std::vector<double> &u_prev) {
  // initial value of the flow
  // TODO: make this more generic for other initial flow values
  const double Q_init = 0;

  for (const auto &v_id : d_graph->get_active_and_connected_vertex_ids(mpi::rank(d_comm))) {
    const auto vertex = d_graph->get_vertex(v_id);

    // exterior boundary
    if (vertex->is_leaf()) {
      const auto edge = d_graph->get_edge(vertex->get_edge_neighbors()[0]);

      const auto &param = edge->get_physical_data();

      // does the vessel point towards the vertex?
      const bool in = edge->is_pointing_to(vertex->get_id());

      const double Q = d_Q_boundary_evaluator(*vertex, *edge);
      const double A = d_A_boundary_evaluator(*vertex, *edge);

      // inflow boundary
      if (vertex->is_inflow_with_fixed_flow()) {
        const double Q_star = (in ? -1 : +1) * vertex->get_inflow_value(t);
        const double A_up = nonlinear::inflow::get_upwinded_A_from_Q(Q, A, in, Q_star, param);

        if (in) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_star;
          d_A_macro_edge_flux_r[edge->get_id()] = A_up;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_star;
          d_A_macro_edge_flux_l[edge->get_id()] = A_up;
        }
      } else if (vertex->is_inflow_with_fixed_pressure()) {
        const double p_up = vertex->get_inflow_value(t);
        const double A_up = nonlinear::get_A_from_p(p_up, param.G0, param.A0);
        const double Q_up = nonlinear::inflow::get_upwinded_Q_from_A(Q, A, (in ? +1 : -1), A_up, param);

        if (in) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_r[edge->get_id()] = A_up;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_l[edge->get_id()] = A_up;
        }
      } else if (vertex->is_free_outflow()) {
        // TODO: make this more generic for other initial flow values
        const double A_init = param.A0;

        double W1, W2;

        if (in) {
          W1 = nonlinear::get_w1_from_QA(Q_init, A_init, param);
          W2 = nonlinear::get_w2_from_QA(Q, A, param);
        } else {
          W1 = nonlinear::get_w1_from_QA(Q, A, param);
          W2 = nonlinear::get_w2_from_QA(Q_init, A_init, param);
        }

        double Q_up = 0, A_up = 0;
        solve_W12(Q_up, A_up, W1, W2, param.G0, param.rho, param.A0);

        if (in) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_r[edge->get_id()] = A_up;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_up;
          d_A_macro_edge_flux_l[edge->get_id()] = A_up;
        }
      } else if (vertex->is_windkessel_outflow() || vertex->is_vessel_tree_outflow()) {
        const bool is_pointing_to = edge->is_pointing_to(vertex->get_id());

        const auto &vertex_dof_map = d_dof_map->get_local_dof_map(*vertex);

        // check that we have dofs assigned for q and p
        assert(vertex_dof_map.num_local_dof() == 1);

        // get Q and p add vessel tip
        const auto &vertex_dofs = vertex_dof_map.dof_indices();

        const auto p_c = u_prev[vertex_dofs[0]];

        const double R1 = calculate_R1(param);

        double Q_out = 0;
        double A_out = 0;
        calculate_windkessel_upwind_values(param, is_pointing_to, R1, Q, A, p_c, Q_out, A_out);

        if (is_pointing_to) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_out;
          d_A_macro_edge_flux_r[edge->get_id()] = A_out;
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_out;
          d_A_macro_edge_flux_l[edge->get_id()] = A_out;
        }
      } else if (vertex->is_nonlinear_characteristic_inflow()) {

        auto &data = vertex->get_nonlinear_characteristic_data();

        double A_r = nonlinear::get_A_from_p(data, data.p);

        std::vector<double> Q_list = {Q, data.q};
        std::vector<double> A_list = {A, A_r};

        std::vector<double> Q_up_list = Q_list;
        std::vector<double> A_up_list = A_list;

        std::vector<VesselParameters> param_list = {
          {param.G0, param.A0, param.rho},
          {data.G0, data.A0, data.rho}};

        std::vector<bool> points_to_vertex_list = {
          edge->is_pointing_to(v_id),
          true};

        solve_at_nfurcation(Q_list, A_list, param_list, points_to_vertex_list, Q_up_list, A_up_list);

        if (edge->is_pointing_to(v_id)) {
          d_Q_macro_edge_flux_r[edge->get_id()] = Q_up_list[0];
          d_A_macro_edge_flux_r[edge->get_id()] = A_up_list[0];
        } else {
          d_Q_macro_edge_flux_l[edge->get_id()] = Q_up_list[0];
          d_A_macro_edge_flux_l[edge->get_id()] = A_up_list[0];
        }
      } else {
        throw std::runtime_error("undefined boundary type!");
      }
    }
  }
}

} // namespace macrocirculation