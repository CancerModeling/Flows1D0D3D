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
  Network(util::BaseModel *model)
      : d_has_network_changed(false), d_model_p(model), d_update_number(0), d_coupled_solver(false), d_comm_p(model->get_comm()), d_procRank(0), d_procSize(0) {

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

  const util::unet::ListStructure<util::unet::VGNode> &get_mesh() const { return VGM; }

  util::unet::ListStructure<util::unet::VGNode> &get_mesh() { return VGM; }

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

  void assemble3D1DSystemForPressure(BaseAssembly
                                       &pres_sys,
                                     BaseAssembly &tum_sys);

  void assemble3D1DSystemForNutrients(BaseAssembly &nut_sys, BaseAssembly &tum_sys);

  void solve3D1DFlowProblem(BaseAssembly
                              &pres_sys,
                            BaseAssembly &tum_sys);

  void solve3D1DNutrientProblem(BaseAssembly &nut_sys, BaseAssembly &tum_sys);

  /*! @brief Solve 1-d system */
  void assembleVGMSystemForPressure(BaseAssembly &pres_sys);

  void solveVGMforPressure(BaseAssembly &pres_sys);

  void assembleVGMSystemForNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  void solveVGMforNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  /** @}*/

  /**
   * @name Update
   */
  /**@{*/

  /*! @brief Update network */
  void updateNetwork(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys);

  void markApicalGrowth();

  void processApicalGrowth();

  void createASingleNode(std::vector<double> new_point, double radius, std::shared_ptr<VGNode> &pointer);

  void createALinkingNode(std::vector<double> new_point, double radius, std::shared_ptr<VGNode> &pointer);

  void linkTerminalVessels();

  bool linkToNearestNetworkNode(std::shared_ptr<VGNode> &pointer);

  bool testCollision(std::vector<double> point);

  bool testIntersection(std::vector<double> point_1, std::vector<double> point_2, double radius, std::shared_ptr<VGNode> &pointer_test);

  void removeRedundantTerminalVessels();

  void markSproutingGrowth();

  void processSproutingGrowth();

  std::vector<double> findNearNetworkNode(std::vector<double> coord, std::vector<double> normal_plane);

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

  /*! @brief Current vessel pressure */
  std::vector<double> P_v;

  /*! @brief Current vessel nutrient */
  std::vector<double> C_v;

  /*! @brief Old vessel nutrient */
  std::vector<double> C_v_old;

  /*! @brief Current 3D1D pressure */
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

  unsigned int d_update_number;
  unsigned int d_update_interval;

  bool d_coupled_solver;

  std::ostringstream oss;

  util::DistributionSample<LogNormalDistribution> d_logNormalDist;
  util::DistributionSample<NormalDistribution> d_normalDist;
  util::DistributionSample<UniformDistribution> d_uniformDist;
};

} // namespace unet

} // namespace util

#endif // UTIL_UNET_NETWORK_H
