////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_LINEARIZED_HEART_TO_BREAST_1D_SOLVER_HPP
#define TUMORMODELS_LINEARIZED_HEART_TO_BREAST_1D_SOLVER_HPP

#include <vector>
#include <string>
#include <memory>
#include <map>

#include <mpi.h>

namespace macrocirculation {

class VesselTipCurrentCouplingData;
class GraphStorage;
class DofMap;
class ImplicitLinearFlowSolver;
class GraphCSVWriter;
class CSVVesselTipWriter;
class GraphPVDWriter;
class Point;
class LinearizedFlowUpwindEvaluator;
class UpwindProviderLinearizedFlow;
class ImplicitTransportSolver;

class LinearizedHeartToBreast1DSolver {
public:
  explicit LinearizedHeartToBreast1DSolver(MPI_Comm comm = MPI_COMM_WORLD);

  void set_path_inflow_pressures(const std::string &path);

  void set_path_geometry(const std::string &path);

  void setup(size_t degree, double tau);

  void solve_flow(double tau, double t);

  void solve_transport(double tau, double t);

  void apply_slope_limiter_transport(double t);

  void write_output(double t);

  void set_output_folder(std::string output_dir);

  std::vector<VesselTipCurrentCouplingData> get_vessel_tip_pressures();

  void update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets);

  void update_vessel_tip_concentrations(const std::map<size_t, double> &concentrations_at_outlets);

private:
  MPI_Comm d_comm;

  size_t d_degree;

  bool d_is_setup;

  std::shared_ptr<GraphStorage> graph;

  std::shared_ptr<DofMap> d_dof_map;

  double d_tau_flow;

  std::string path_linear_geometry{"data/1d-meshes/coarse-breast-geometry-with-extension.json"};

  std::string path_inflow_pressures{"data/1d-input-pressures/from-33-vessels-with-small-extension.json"};

  std::string output_folder_name{"output"};
  std::string filename_csv{"linearized_heart_to_breast_1d_solution"};
  std::string filename_csv_tips{"linearized_heart_to_breast_1d_solution_tips"};
  std::string filename_pvd{"linearized_heart_to_breast_1d_solution"};

  std::shared_ptr<ImplicitLinearFlowSolver> solver;

  std::shared_ptr<GraphCSVWriter> csv_writer;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer;
  std::shared_ptr<GraphPVDWriter> graph_pvd_writer;

  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;
  std::vector<double> c_vertex_values;

  std::vector<double> vessel_ids;
  std::vector<double> vessel_radii;
  std::vector<double> vessel_A;

  std::shared_ptr<LinearizedFlowUpwindEvaluator> upwind_evaluator;
  std::shared_ptr<UpwindProviderLinearizedFlow> upwind_provider;
  std::shared_ptr<ImplicitTransportSolver> transport_solver;

  void setup_graphs();

  void setup_solver(size_t degree, double tau);

  void setup_solver_flow(size_t degree, double tau);

  void setup_solver_transport(size_t degree);

  void setup_output();
};

}

#endif //TUMORMODELS_LINEARIZED_HEART_TO_BREAST_1D_SOLVER_HPP
