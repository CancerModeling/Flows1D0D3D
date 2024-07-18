////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2022 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "dof_map.hpp"
#include "graph_partitioner.hpp"
#include "graph_pvd_writer.hpp"
#include "graph_csv_writer.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "simple_nonlinear_solver.hpp"
#include "vessel_formulas.hpp"
#include "embedded_graph_reader.hpp"
#include "petsc.h"
#include "communication/mpi.hpp"
#include "interpolate_to_vertices.hpp"
#include "quantities_of_interest.hpp"

namespace macrocirculation {

SimpleNonlinearSolver::SimpleNonlinearSolver(MPI_Comm comm, const std::string& filepath_mesh, const std::string& filepath_output, const std::string& name, double dt)
// TODO: Manage parallelism better
  : comm(comm)
  , tau(dt)
  , t(0)
  , graph(std::make_shared<GraphStorage>())
  , pvd_writer(std::make_shared<GraphPVDWriter>(comm, filepath_output, name))
  , csv_writer(std::make_shared<GraphCSVWriter>(comm, filepath_output, name, graph)) 
  , num_coupling_points(0)
{
  EmbeddedGraphReader graph_reader;
  graph_reader.append(filepath_mesh, *graph);

  for (int i = 0; i <= 100; i+=1) // more than 100 coupling vertices are not realistic
  {
    if (graph->has_named_vertex("coupling_" + std::to_string(i) + "_outer"))
    {
      num_coupling_points += 1;

      v_coupling_outer.push_back(graph->find_vertex_by_name("coupling_" + std::to_string(i) + "_outer"));
      v_coupling_inner.push_back(graph->find_vertex_by_name("coupling_" + std::to_string(i) + "_inner"));
      edge.push_back(graph->find_edge_by_name("vessel_coupling_" + std::to_string(i)));

      if ( v_coupling_outer.back()->is_leaf() )
      {
        const double C = linear::get_C(edge.back()->get_physical_data());
        const double L = linear::get_L(edge.back()->get_physical_data());
        const double G0 = edge.back()->get_physical_data().G0;
        const double A0 = edge.back()->get_physical_data().A0;
        const double rho = edge.back()->get_physical_data().rho;
        const bool ptv = edge[0]->is_pointing_to( v_coupling_outer.back()->get_id() );
        v_coupling_outer.back()->set_to_nonlinear_characteristic_inflow(G0, A0, rho, true, 0, 0);
      }
      else
      {
        throw std::runtime_error("Outer coupling point must be a leaf node.");
      }
    }
    else
    {
      break;
    }
  }

  if (graph->has_named_vertex("Inflow"))
  {
    auto v0 = graph->find_vertex_by_name("Inflow");
    v0->set_to_inflow_with_fixed_flow([](double t) { return 2 * std::abs(std::sin(M_PI * t)); });
  }

  /*
  if (graph->has_named_vertex("outflow"))
  {
    auto v3 = graph->find_vertex_by_name("outflow");
    v3->set_to_free_outflow();
  }
  */
}

SimpleNonlinearSolver::SimpleNonlinearSolver(const std::string& filepath_mesh, const std::string& filepath_output, const std::string& name, double dt)
// TODO: Manage parallelism better
  : SimpleNonlinearSolver(PETSC_COMM_WORLD, filepath_mesh, filepath_output, name, dt) {};

void SimpleNonlinearSolver::set_tau(double ptau) {
  tau = ptau;
}

void SimpleNonlinearSolver::set_inflow(const std::function< double(double)> & fun) {
  if (graph->has_named_vertex("Inflow"))
  {
    auto v0 = graph->find_vertex_by_name("Inflow");
    v0->set_to_inflow_with_fixed_flow(fun);
  }
}

void SimpleNonlinearSolver::set_outflow_rcr(const double r, const double c) {
  if (graph->has_named_vertex("Outflow"))
  {
    auto v3 = graph->find_vertex_by_name("Outflow");
    v3->set_to_windkessel_outflow(r, c);
  }
}

void SimpleNonlinearSolver::setup(){
  const size_t degree = 0;

  graph->finalize_bcs();

  naive_mesh_partitioner(*graph, comm);

  dof_map = std::make_shared<DofMap>(graph->num_vertices(), graph->num_edges());
  dof_map->create(comm, *graph, 2, degree, false);

  solver = std::make_shared<ExplicitNonlinearFlowSolver>(comm, graph, dof_map, degree);
  solver->use_explicit_euler_method();

  csv_writer->add_setup_data(dof_map, solver->A_component, "A");
  csv_writer->add_setup_data(dof_map, solver->Q_component, "Q");
  csv_writer->setup();
}

void SimpleNonlinearSolver::solve() {
  if (!graph->bcs_finalized())
    setup();

  t += tau;
  solver->solve(tau, t);
}

double SimpleNonlinearSolver::get_current_t() const { return t; }

SimpleNonlinearSolver::Result SimpleNonlinearSolver::get_result(const Vertex &vertex, const Edge &edge) {
  double p, q, a;
  if (edge.rank() == mpi::rank(comm))
  {
    solver->get_1d_pq_values_at_vertex(vertex, edge, p, q);
    const auto &phys_data = edge.get_physical_data();
    a = nonlinear::get_A_from_p(p, phys_data.G0, phys_data.A0);
    p *= 1e3; // mixed units to cgs

    // vertex is the inner vertex. if the edge points to it, it points to the wrong direction
    if (edge.is_pointing_to(vertex.get_id()))
      q *= -1;
  }
  // communicate everywhere:
  double data[3]  = {p, q, a};
  MPI_Bcast(data, 3, MPI_DOUBLE, edge.rank(), comm);
  p = data[0]; q = data[1]; a = data[2];
  // return
  return {a, p, q};
}

void SimpleNonlinearSolver::set_result(int outlet, double p, double q) {
  p /= 1e3; // cgs to mixed units
  if (outlet >= num_coupling_points)
    throw std::runtime_error("unknown vertex value");
  v_coupling_outer[outlet]->update_nonlinear_characteristic_inflow(p, q);
}

SimpleNonlinearSolver::Result SimpleNonlinearSolver::get_result(int outlet) {
  if (outlet >= num_coupling_points)
    throw std::runtime_error("unknown vertex value");
  return get_result(*v_coupling_inner[outlet], *edge[outlet]);
}

std::vector<std::array<double, 3>> SimpleNonlinearSolver::get_points()
{
  // TODO: Generalize this for arbitrary orientations
  const auto v1 = edge[0]->get_embedding_data().points.front();
  const auto v2 = edge[0]->get_embedding_data().points.back();
  const auto v3 = edge[1]->get_embedding_data().points.front();
  const auto v4 = edge[1]->get_embedding_data().points.back();
  return {
    {v1.x, v1.y, v1.z},
    {v2.x, v2.y, v2.z},
    {v3.x, v3.y, v3.z},
    {v4.x, v4.y, v4.z},
  };
}

SimpleNonlinearSolver::Result SimpleNonlinearSolver::get_result_outer(int outlet) {
  if (outlet >= num_coupling_points)
    throw std::runtime_error("unknown vertex value");

  return get_result(*v_coupling_outer[outlet], *edge[outlet]);
}

void SimpleNonlinearSolver::write() {
  std::vector<Point> points;
  std::vector<double> a_vertex_values;
  std::vector<double> q_vertex_values;
  std::vector<double> p_vertex_values;
  std::vector<double> vertex_ids;
  std::vector<double> vessel_ids;
  std::vector<double> A0;
  std::vector<double> r_vertex_values;
  interpolate_to_vertices(comm, *graph, *dof_map, solver->Q_component, solver->get_solution(), points, q_vertex_values);
  interpolate_to_vertices(comm, *graph, *dof_map, solver->A_component, solver->get_solution(), points, a_vertex_values);
  fill_with_vessel_id(comm, *graph, points, vessel_ids);
  fill_with_vertex_id(MPI_COMM_WORLD, *graph, points, vertex_ids);
  fill_with_vessel_A0(comm, *graph, points, A0);
  calculate_nonlinear_p(comm, *graph, a_vertex_values, points, p_vertex_values);

  pvd_writer->set_points(points);
  pvd_writer->add_vertex_data("a", a_vertex_values);
  pvd_writer->add_vertex_data("q", q_vertex_values);
  pvd_writer->add_vertex_data("p", p_vertex_values);
  pvd_writer->add_vertex_data("vessel_ids", vessel_ids);
  pvd_writer->add_vertex_data("vertex_ids", vertex_ids);
  pvd_writer->add_vertex_data("a0", A0);
  pvd_writer->write(t);

  csv_writer->add_data("A", solver->get_solution());
  csv_writer->add_data("Q", solver->get_solution());
  csv_writer->write(t);
}

size_t SimpleNonlinearSolver::get_num_coupling_points() const { return num_coupling_points; }

} // namespace macrocirculation