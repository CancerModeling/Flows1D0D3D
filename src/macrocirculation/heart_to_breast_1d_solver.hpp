////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_HEART_TO_BREAST_1_D_SOLVER_HPP
#define TUMORMODELS_HEART_TO_BREAST_1_D_SOLVER_HPP

#include "mpi.h"
#include <memory>

// include for points
// TODO: move points into its own file
#include "graph_storage.hpp"

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class NonlinearLinearCoupling;
class ImplicitLinearFlowSolver;
class ExplicitNonlinearFlowSolver;
class CoupledExplicitImplicit1DSolver;
class GraphCSVWriter;
class GraphPVDWriter;
class CSVVesselTipWriter;
class TipVertexDofIntegrator;

struct VesselTipCouplingData {
  Point p;
  double p_art;
  double p_ven;
  double R2_art;
  double R2_cap;
};

class HeartToBreast1DSolver {
public:
  explicit HeartToBreast1DSolver(MPI_Comm comm);

  void setup(size_t degree, double tau);

  void solve();

  void write_output();

  double get_time() const;

  void start_0d_pressure_integrator();

  std::vector<VesselTipCouplingData> stop_0d_pressure_integrator();

private:
  MPI_Comm d_comm;

  size_t d_degree;

  std::shared_ptr<GraphStorage> graph_nl;
  std::shared_ptr<GraphStorage> graph_li;

  std::shared_ptr<NonlinearLinearCoupling> coupling;

  double d_tau;

  double t;

  std::string path_nonlinear_geometry{"data/meshes/network-33-vessels-extended.json"};
  std::string path_linear_geometry{"data/meshes/coarse-network-geometry.json"};

  std::string path_boundary_nonlinear{"data/meshes/boundary-combined-geometry-nonlinear-part.json"};
  std::string path_boundary_linear{"data/meshes/boundary-combined-geometry-linear-part.json"};

  std::string output_folder_name{"output"};
  std::string filename_csv_nl{"combined_geometry_solution_nl"};
  std::string filename_csv_li{"combined_geometry_solution_li"};
  std::string filename_csv_tips{"combined_geometry_solution"};
  std::string filename_pvd{"combined_geometry_solution_tips"};

  std::shared_ptr<CoupledExplicitImplicit1DSolver> solver;

  std::shared_ptr<TipVertexDofIntegrator> integrator;

  size_t last_arterial_tip_index{3};
  size_t capillary_tip_index{4};
  size_t first_vene_tip_index{5};

  std::shared_ptr<GraphCSVWriter> csv_writer_nl;
  std::shared_ptr<GraphCSVWriter> csv_writer_li;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer;
  std::shared_ptr<GraphPVDWriter> graph_pvd_writer;

  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;

  std::vector<double> vessel_ids_li;

  bool d_integrator_running;

  ImplicitLinearFlowSolver &get_solver_li();
  ExplicitNonlinearFlowSolver &get_solver_nl();

  void setup_graphs();

  void setup_solver(size_t degree, double tau);

  void setup_output();
};

} // namespace macrocirculation


#endif //TUMORMODELS_HEART_TO_BREAST_1_D_SOLVER_HPP
