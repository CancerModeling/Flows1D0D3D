////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNET_NETWORK_H
#define UTIL_UNET_NETWORK_H

// Libmesh
#include "random_dist.hpp"
#include "utils.hpp"

// gmm dependencies
#include "gmm.h"

// custom data structures
#include "list_structure.hpp"
#include "nodes.hpp"

// network utility
#include "netUtil.hpp"

// model class
#include "umodel/model.hpp"

// assembly class
#include "usystem/abstraction.hpp"

namespace util {

inline double pressure_initial_guess_constant(double p) { return p; }
inline double pressure_initial_guess_95_percent(double p) { return 0.95 * p; }

/*!
 * @brief Namespace for 1D network
 */
namespace unet {

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Network {

public:
  /*! @brief Constructor */
  explicit Network(util::BaseModel *model)
      : d_has_network_changed(false),
        d_model_p(model),
        d_update_number(0),
        d_coupled_solver(false),
        d_extrapolate_nutrients_at_tips(false),
        d_comm_p(model->get_comm()),
        d_procRank(0),
        d_procSize(0),
        total_added_length(0.0),
        total_removed_length(0.0),
        total_added_volume(0.0),
        total_removed_volume(0.0),
        d_pressure_boundary_initial_guess(model->get_input_deck().d_pressure_initial_guess_95_percent ? pressure_initial_guess_95_percent : pressure_initial_guess_constant) {

    // initialize random distribution samplers
    const auto &input = d_model_p->get_input_deck();

    if (d_comm_p) {
      d_procSize = d_comm_p->size();
      d_procRank = d_comm_p->rank();
    }

    // TAG: Random
    d_logNormalDist.init(input.d_log_normal_mean, input.d_log_normal_std_dev, input.d_seed);
    d_normalDist.init(0., 1., input.d_seed);
    d_uniformDist.init(0., 1., input.d_seed);
  }

  const ListStructure<VGNode> &get_mesh() const { return VGM; }

  ListStructure<VGNode> &get_mesh() { return VGM; }

  /**
   * @name Vector-Copy-Functions
   * @brief Utility functions to copy data into and out of the network.
   */
  /**@{*/

  /*! @brief Unmarks all the nodes in the network by setting the node_marked attribute to false. */
  void unmark_network_nodes();

  /*! @brief Marks all the nodes connected to an initial node by setting the node_marked attribute to true. */
  void mark_nodes_connected_with_initial_nodes();

  /*! @brief Adds the volume and length of all vessels in contact with unmarked vertices.
   *         This function call normally precedes a call which deletes the unmarked to gather statistics about edge removal.
   */
  void add_lengths_and_volumes_of_unmarked_network(double &length, double &volume);

  /*! @brief Gets the overall volume and length of all vessels in the network. */
  void get_length_and_volume_of_network(double &length, double &volume);

  /*! @brief Returns the number of bifurcation in the network. */
  int get_number_of_bifurcations();

  /*! @brief Deletes nodes with a node_marked attribute set to true. */
  void delete_unmarked_nodes();

  /*! @brief Reenumerates the network dofs. Must be called after deleting nodes. */
  void reenumerate_dofs();

  /*! @brief Removes all the nodes which are not connected with the initial network. */
  void delete_unconnected_nodes();

  /*! @brief Removes all the nodes which were created for sprouting, but whose vessels were removed in the process. */
  void delete_old_sprouters();

  /*! @brief Returns the maximum pressure value in 3D. */
  double get_maximum_pressure_value_3D() const;

  /*! @brief Sets the boundary condition of all vessel tips not present in the original mesh to Neumann. */
  void set_bc_of_added_vessels_to_neumann();

  /*! @brief Sets the edge_touched and sprouting_edge flags to false. */
  void reset_edge_flags();

  /*! @brief Resizes the matrices of the direct solver, after the network has changed. */
  void resize_matrices_direct_solver();

  /*! @brief Resizes the matrices of the coupled solver, after the network has changec. */
  void resize_matrices_coupled_solver();

  /*! @brief Copies the dof values stored on the network nodes (pressure and nutrients) into the vectors. */
  void copy_network_to_vectors();

  /**
   * @name Input-output
   */
  /**@{*/

  /*! @brief Create mesh for 1d network */
  void create_initial_network();

  void readData(std::vector<std::vector<double>> &vertices,
                std::vector<double> &pressures, std::vector<double> &radii,
                std::vector<std::vector<unsigned int>> &elements);

  void transferDataToVGM(std::vector<std::vector<double>> &vertices,
                         std::vector<double> &pressures,
                         std::vector<double> &radii,
                         std::vector<std::vector<unsigned int>> &elements);

  void printDataVGM();

  void writeDataToVTK_3D(std::vector<double> P_3D, int N_3D, double h_3D);

  void writeDataToVTK3D_Pressure(std::vector<double> P_3D, std::vector<std::vector<double>> V_3D, int N_3D, double h_3D, int timeStep);

  void writeDataToVTK3D_Nutrients(std::vector<double> C_3D, int N_3D, double h_3D, int timeStep);

  /** @}*/

  /**
   * @name Solve
   */
  /**@{*/

  /*! @brief Assembles the fully coupled 1D-3D pressure system */
  void assemble3D1DSystemForPressure(BaseAssembly &pres_sys, BaseAssembly &tum_sys);

  /*! @brief Assembles the fully coupled 1D-3D nutrient system */
  void assemble3D1DSystemForNutrients(BaseAssembly &nut_sys, BaseAssembly &tum_sys);

  /*! @brief Solves the fully coupled 1D-3D pressure problem */
  void solve3D1DFlowProblem(BaseAssembly &pres_sys, BaseAssembly &tum_sys);

  /*! @brief Solves the fully coupled 1D-3D nutrient problem */
  void solve3D1DNutrientProblem(BaseAssembly &nut_sys, BaseAssembly &tum_sys);

  /*! @brief Assembles the 1D pressure system */
  void assembleVGMSystemForPressure(BaseAssembly &pres_sys);

  /*! @brief Solves one iteration of the decoupled 1D-3D pressure problem */
  void solveVGMforPressure(BaseAssembly &pres_sys);

  /*! @brief Assembles the 1D nutrient system */
  void assembleVGMSystemForNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  /*! @brief Solves one iteration of the decoupled 1D-3D nutrient problem */
  void solveVGMforNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  /** @}*/

  /**
   * @name Update
   */
  /**@{*/

  /*! @brief Update network */
  void updateNetwork(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys);

  /*! @brief Marks all candidates for apical growth by setting the apicalGrowth flag to true. */
  void markApicalGrowth();

  void processApicalGrowth();

  void createASingleNode(std::vector<double> new_point, double radius, std::shared_ptr<VGNode> &pointer);

  /*! @brief Links the terminal vessels with other nearby vessels. */
  void linkTerminalVessels();

  bool linkToNearestNetworkNode(std::shared_ptr<VGNode> &pointer);

  bool testCollision(std::vector<double> point);

  bool testIntersection(const std::vector<double>& point_1, const std::vector<double>& point_2, double radius, std::shared_ptr<VGNode> &pointer_test);

  void removeRedundantTerminalVessels();

  /*! @brief Marks all edge candidates for sprouting growth.
   *
   */
  void markSproutingGrowth();

  void processSproutingGrowth();

  std::vector<double> findNearNetworkNode(std::vector<double> coord, std::vector<double> normal_plane);

  /*! @brief Adapts the radii of all vessels. */
  void adaptRadius();

  /** @}*/

  /**
   * @name Utility functions
   */
  /**@{*/

  void update_old_concentration() {
    if (d_coupled_solver)
      phi_sigma_old = phi_sigma;
    else
      C_v_old = C_v;
  };

  void update_old_concentration_3D1D() { phi_sigma_old = phi_sigma; };

  void refine1DMesh();

  double getDirichletValue(std::vector<double> center_face, double L_p, double radius);

  double getK1D(double s, double L_p, double radius);

  void rescaleSecombData(std::vector<std::vector<double>> &vertices, std::vector<double> &pressures, std::vector<double> &radii, double epsilon);

  double sourceTermTAFTwoVessels(std::vector<double> coord);

  std::vector<util::unet::ElemWeights>
  compute_elem_weights_at_node(std::shared_ptr<util::unet::VGNode> &pointer)
    const;

  /*! @brief Compute intersecting element and weight in tumor domain */
  void compute_elem_weights();

  /*! @brief Compute qoi */
  std::vector<double> compute_qoi();

  unsigned int get_assembly_cases_pres(const std::shared_ptr<VGNode> &pointer, const double &identify_vein_pres) const;
  unsigned int get_assembly_cases_nut(const std::shared_ptr<VGNode> &pointer, const double &identify_vein_pres) const;

  std::string get_assembly_cases_pres_str(const std::shared_ptr<VGNode> &pointer, const double &identify_vein_pres) const;
  std::string get_assembly_cases_nut_str(const std::shared_ptr<VGNode> &pointer, const double &identify_vein_pres) const;

  void prepare_and_communicate_network();
  void update_and_communicate_bdry_flag();
  void check_vessel_length();

  void removeSingleEdges();

  /** @}*/

  /*! @brief Did we change network from last step */
  bool d_has_network_changed;

  /*! @brief Reference to model class */
  util::BaseModel *d_model_p;

  /*! @brief Pointer to communicator */
  Parallel::Communicator *d_comm_p;

  /*! @brief Rank of processor */
  unsigned int d_procRank;

  /*! @brief Size of processors */
  unsigned int d_procSize;

  /*! @brief 1-D mesh */
  util::unet::ListStructure<util::unet::VGNode> VGM;

  /*! @brief number of vertices in 1D network */
  unsigned int d_numVertices;

  /*! @brief number of segments in 1D network */
  unsigned int d_numSegments;

  /*! @brief Vertices data
   *
   * (x-coord, y-coord, z-coord) and boundry flag
   *
   */
  std::vector<double> d_vertices;
  std::vector<unsigned int> d_vertexBdFlag;

  /*! @brief Segments data
   *
   * (node 1, node 2) and radius
   *
   */
  std::vector<unsigned int> d_segments;
  std::vector<double> d_segmentData;


  /*! @brief System matrix for vessel pressure */
  gmm::row_matrix<gmm::wsvector<double>> A_VGM;

  /*! @brief System matrix for vessel nutrient */
  //gmm::row_matrix<gmm::wsvector<double>> Ac_VGM;

  /*! @brief System matrix for 3D1D pressure flow model*/
  gmm::row_matrix<gmm::wsvector<double>> A_flow_3D1D;

  /*! @brief System matrix for 3D1D nutrient model*/
  gmm::row_matrix<gmm::wsvector<double>> A_nut_3D1D;

  /*! @brief System force for vessel pressure */
  std::vector<double> b;

  /*! @brief System force for vessel nutrient */
  std::vector<double> b_c;

  /*! @brief Right hand side pressure 3D1D */
  std::vector<double> b_flow_3D1D;

  /*! @brief Right hand side nutrients 3D1D */
  std::vector<double> b_nut_3D1D;

  /*! @brief Current vessel pressures used by the _direct_ solver. */
  std::vector<double> P_v;

  /*! @brief Current vessel nutrients used by the _direct_ solver. */
  std::vector<double> C_v;

  /*! @brief Old vessel nutrients used by the _direct_ solver. */
  std::vector<double> C_v_old;

  /*! @brief Current pressure values, used by the _coupled_ solver. */
  std::vector<double> P_3D1D;

  /*! @brief Current 3D pressure */
  std::vector<double> P_3D;

  /*! @brief Current 3D1D nutrient concentration */
  std::vector<double> phi_sigma;

  /*! @brief Old 3D1D nutrient concentration */
  std::vector<double> phi_sigma_old;

  /*! @brief Current 3D nutrient concentration */
  std::vector<double> phi_sigma_3D;

  /*! @brief Current TAF concentration */
  std::vector<double> phi_TAF;

  double mu;

  double D_v;

  double D_v_3D;

  double D_TAF;

  double osmotic_sigma;

  std::string scenario;

  double h_3D;

  int N_3D;

  int N_tot_3D;

  double K_3D;
  double L_x;

  double L_p;
  double L_s;

  /*! @brief Counts how much vessel length is added during a growth step. */
  double total_added_length;

  /*! @brief Counts how much vessel length is removed during a growth step. */
  double total_removed_length;

  /*! @brief Counts how much vessel vessel volume is added during a growth step. */
  double total_added_volume;

  /*! @brief Counts how much vessel vessel volume is removed during a growth step. */
  double total_removed_volume;

  unsigned int d_update_number;
  unsigned int d_update_interval;

  bool d_coupled_solver;

  bool d_extrapolate_nutrients_at_tips;

  std::ostringstream oss;

  util::DistributionSample<LogNormalDistribution> d_logNormalDist;
  util::DistributionSample<NormalDistribution> d_normalDist;
  util::DistributionSample<UniformDistribution> d_uniformDist;

  /*! Provides an initial guess for the boundary value at the new grown vessel tips. */
  std::function<double(double)> d_pressure_boundary_initial_guess;
};

} // namespace unet

} // namespace util

#endif // UTIL_UNET_NETWORK_H
