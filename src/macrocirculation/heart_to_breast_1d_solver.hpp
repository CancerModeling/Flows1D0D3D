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

struct VesselTipAverageCouplingData {
  /*! @brief Coordinates of the vessel tip. */
  Point p;

  /*! @brief Vertex id at the vessel tip. */
  size_t vertex_id;

  /*! @brief Pressure value in the arterioles in [Ba]. */
  double p_art;

  /*! @brief Pressure value in the venules in [Ba]. */
  double p_ven;

  /*! @brief Second resistance of the arteriole compartment in [cm^{-4} g s^{-1}]. */
  double R2_art;

  /*! @brief Second resistance of the capillary compartment in [cm^{-4} g s^{-1}]. */
  double R2_cap;
};

struct VesselTipCurrentCouplingData {
  /*! @brief Coordinates of the vessel tip. */
  Point p;

  /*! @brief Vertex id at the vessel tip. */
  size_t vertex_id;

  /*! @brief Pressure value in the arterioles in [Ba]. */
  double pressure;

  /*! @brief Second resistance of the arteriole compartment in [cm^{-4} g s^{-1}]. */
  double R2;
};

/*! @brief Denotes Different boundary conditions. */
enum class BoundaryModel {
  /*! @brief RCR-models, which are chained one after the other. */
  DiscreteRCRChain,
  /*! @brief RCR-models, which form a tree like structure. */
  DiscreteRCRTree
};

class HeartToBreast1DSolver {
public:
  explicit HeartToBreast1DSolver(MPI_Comm comm);

  void setup(size_t degree, double tau, BoundaryModel boundary_model);

  void solve();

  void write_output();

  double get_time() const;

  void start_0d_pressure_integrator();

  /*! @brief Returns the *averaged* pressures at the vessel tips in CGS units. */
  std::vector<VesselTipAverageCouplingData> stop_0d_pressure_integrator();

  /*! @brief Returns the *current* pressures at the vessel tips in CGS units. */
  std::vector<VesselTipCurrentCouplingData> get_vessel_tip_pressures();

  /*! @brief Updates the boundary pressures at the vessel outlets in the linearized regime.
   *         The keys of the given map are the vertex ids of the boundary values, while the values correspond to the pressure values.
   *         On one processor all the keys for the owned outlets have to be present in the map, otherwise an exception is thrown.
   */
  void update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets);

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
  std::string filename_csv_tips{"combined_geometry_solution_tips"};
  std::string filename_pvd{"combined_geometry_solution"};

  std::shared_ptr<CoupledExplicitImplicit1DSolver> solver;

  std::shared_ptr<TipVertexDofIntegrator> integrator;

  size_t last_arterial_tip_index{7};
  size_t capillary_tip_index{8};
  size_t first_vene_tip_index{9};

  std::shared_ptr<GraphCSVWriter> csv_writer_nl;
  std::shared_ptr<GraphCSVWriter> csv_writer_li;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer;
  std::shared_ptr<GraphPVDWriter> graph_pvd_writer;

  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;

  std::vector<double> vessel_ids_li;
  std::vector<double> vessel_radii_li;

  bool d_integrator_running;

  ImplicitLinearFlowSolver &get_solver_li();
  ExplicitNonlinearFlowSolver &get_solver_nl();

  void setup_graphs(BoundaryModel boundary_model);

  void setup_solver(size_t degree, double tau);

  void setup_output();
};

} // namespace macrocirculation


#endif //TUMORMODELS_HEART_TO_BREAST_1_D_SOLVER_HPP
