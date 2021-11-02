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
class UpwindProvider;
class ImplicitTransportSolver;
class UpwindProviderLinearizedFlow;
class UpwindProviderNonlinearFlow;
class LinearizedFlowUpwindEvaluator;
class NonlinearFlowUpwindEvaluator;
class DofMap;
class PetscVec;

struct VesselTipCurrentCouplingData {
  /*! @brief Coordinates of the vessel tip. */
  Point p;

  /*! @brief Vertex id at the vessel tip. */
  size_t vertex_id;

  /*! @brief Pressure value in the arterioles in [Ba]. */
  double pressure;

  /*! @brief Concentration value at the tip in [mmol cm^{-3}]. */
  double concentration;

  /*! @brief Second resistance of the arteriole compartment in [cm^{-4} g s^{-1}]. */
  double R2;

  /*! @brief The radius of the last vessel tree segment in [cm]. */
  double radius;

  /*! @brief The number of levels for this tree. */
  size_t level;
};

class HeartToBreast1DSolver {
public:
  explicit HeartToBreast1DSolver(MPI_Comm comm);

  void set_path_inflow_pressures(const std::string& path);

  void set_path_nonlinear_geometry(const std::string& path);

  void set_path_linear_geometry(const std::string& path);

  void set_path_coupling_conditions(const std::string& path);

  void setup(size_t degree, double tau);

  void solve_flow(double tau, double t);

  void solve_transport(double tau, double t);

  void apply_slope_limiter_transport(double t);

  void write_output(double t);

  /*! @brief Returns the *current* pressures at the vessel tips in CGS units. */
  std::vector<VesselTipCurrentCouplingData> get_vessel_tip_pressures();

  /*! @brief Updates the boundary pressures at the vessel outlets in the linearized regime.
   *         The keys of the given map are the vertex ids of the boundary values, while the values correspond to the pressure values.
   *         On one processor all the keys for the owned outlets have to be present in the map, otherwise an exception is thrown.
   *         The pressure unit has to be [Ba].
   */
  void update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets);

  /*! @brief Updates the boundary concentrations at the vessel outlets of the linearized regime.
   *         The keys of the given map are the vertex ids of the boundary values, while the values correspond to the concentration values.
   *         The concentration is given in [mmol cm^{-3}]
   */
  void update_vessel_tip_concentrations(const std::map<size_t, double> &concentrations_at_outlets);

  /*! @brief Set output folder name. */
  void set_output_folder(std::string output_dir);

private:
  MPI_Comm d_comm;

  size_t d_degree;

  bool d_is_setup;

  std::shared_ptr<GraphStorage> graph_nl;
  std::shared_ptr<GraphStorage> graph_li;

  std::shared_ptr<NonlinearLinearCoupling> coupling;

  double d_tau_flow;

  std::string path_nonlinear_geometry{"data/1d-meshes/33-vessels-with-small-extension.json"};
  std::string path_linear_geometry{"data/1d-meshes/coarse-breast-geometry-with-extension.json"};

  std::string path_inflow_pressures{""};

  std::string path_coupling_conditions{"data/1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"};

  std::string output_folder_name{"output"};
  std::string filename_csv_nl{"heart_to_breast_1d_solution_nl"};
  std::string filename_csv_li{"heart_to_breast_1d_solution_li"};
  std::string filename_csv_tips_nl{"heart_to_breast_1d_solution_tips_nl"};
  std::string filename_csv_tips_li{"heart_to_breast_1d_solution_tips_li"};
  std::string filename_pvd{"heart_to_breast_1d_solution"};

  std::shared_ptr<CoupledExplicitImplicit1DSolver> solver;

  std::shared_ptr<TipVertexDofIntegrator> integrator;

  std::shared_ptr<GraphCSVWriter> csv_writer_nl;
  std::shared_ptr<GraphCSVWriter> csv_writer_li;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer_nl;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer_li;
  std::shared_ptr<GraphPVDWriter> graph_pvd_writer;

  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;
  std::vector<double> c_vertex_values;

  std::vector<double> vessel_ids_li;
  std::vector<double> vessel_radii_li;
  std::vector<double> vessel_A_li;

  std::shared_ptr<NonlinearFlowUpwindEvaluator> upwind_evaluator_nl;
  std::shared_ptr<LinearizedFlowUpwindEvaluator> upwind_evaluator_li;

  std::shared_ptr<UpwindProviderNonlinearFlow> upwind_provider_nl;
  std::shared_ptr<UpwindProviderLinearizedFlow> upwind_provider_li;

  std::shared_ptr<ImplicitTransportSolver> transport_solver;

  ImplicitLinearFlowSolver &get_solver_li();
  ExplicitNonlinearFlowSolver &get_solver_nl();

  void setup_graphs();

  void setup_solver(size_t degree, double tau);

  void setup_solver_flow(size_t degree, double tau);

  void setup_solver_transport(size_t degree);

  void setup_output();
};

std::map<size_t, double> get_vessel_tip_dof_values(MPI_Comm comm, const GraphStorage &graph, const DofMap &dof_map, const PetscVec &u);

} // namespace macrocirculation


#endif //TUMORMODELS_HEART_TO_BREAST_1_D_SOLVER_HPP
