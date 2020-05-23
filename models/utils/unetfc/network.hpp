////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNETFC_NETWORK_H
#define UTIL_UNETFC_NETWORK_H

// Libmesh
#include "utils.hpp"

// gmm dependencies
#include "gmm.h"

// custom data structures
#include "list_structure.hpp"
#include "nodes.hpp"

// model class
#include "umodel/model.hpp"

// assembly class
#include "usystem/abstraction.hpp"

namespace util {

/*!
 * @brief Namespace for 1D network
 */
namespace unetfc {

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Network {

public:
  /*! @brief Constructor */
  Network(util::BaseModel *model)
      : d_is_network_changed(false), d_model_p(model), d_update_number(0) {}

  const util::unetfc::ListStructure<util::unetfc::VGNode> &get_mesh() const { return
  VGM; }

  util::unetfc::ListStructure<util::unetfc::VGNode> &get_mesh() { return VGM; }

  /**
   * @name Input-output
   */
  /**@{*/

  /*! @brief Create mesh for 1d network */
  void create_initial_network();

  void readData( std::vector<std::vector<double>> &vertices,
                 std::vector<double> &pressures, std::vector<double> &radii,
                 std::vector<std::vector<unsigned int>> &elements );

  void transferDataToVGM( std::vector<std::vector<double>> &vertices,
                          std::vector<double> &pressures,
                          std::vector<double> &radii,
                          std::vector<std::vector<unsigned int>> &elements );

  void printDataVGM();

  void writeDataToVTKTimeStep_VGM( int timeStep );

  void writeDataToVTK_3D( std::vector<double> P_3D, int N_3D, double h_3D );

  void writeDataToVTK3D_Pressure(std::vector<double> P_3D, std::vector< std::vector<double> > V_3D, int N_3D, double h_3D, int timeStep);

  void writeDataToVTK3D_Nutrients(std::vector<double> C_3D, int N_3D, double h_3D, int timeStep);

  /** @}*/

  /**
   * @name Solve
   */
  /**@{*/

  void assemble3D1DSystemForPressure(BaseAssembly
                                     &pres_sys, BaseAssembly &tum_sys);

  void assemble3D1DSystemForNutrients(BaseAssembly &nut_sys, BaseAssembly &tum_sys);

   void solve3D1DFlowProblem(BaseAssembly
   &pres_sys, BaseAssembly &tum_sys);

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
  void update_network(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys) {}

  void updateNetwork();

  void markApicalGrowth();

  void processApicalGrowth();

  void createASingleNode( std::vector<double> new_point, double radius, std::shared_ptr<VGNode>& pointer );

  void createALinkingNode( std::vector<double> new_point, double radius, std::shared_ptr<VGNode>& pointer );

  void linkTerminalVessels();

  bool testCollision( std::vector<double> point );

  bool testIntersection( std::vector<double> point_1, std::vector<double> point_2, std::vector<double>& new_point_link, double radius );

  void removeRedundantTerminalVessels();

  void markSproutingGrowth();

  void processSproutingGrowth();

  std::vector<double> findNearNetworkNode( std::vector<double> coord, std::vector<double> normal_plane );

  /** @}*/

  /**
   * @name Utility functions
   */
  /**@{*/

  void update_old_concentration() { C_v_old = C_v; };

  void update_old_concentration_3D1D() { phi_sigma_old = phi_sigma; };

  void refine1DMesh();

  double getDirichletValue( std::vector< double > center_face, double L_p, double radius );

  double getK1D( double s, double L_p, double radius );

  void rescaleSecombData( std::vector< std::vector<double> >& vertices, std::vector<double>& pressures, std::vector<double>& radii, double epsilon );

  double sourceTermTAFTwoVessels( std::vector<double> coord );

  std::vector<util::unetfc::ElemWeights>
  compute_elem_weights_at_node(std::shared_ptr<util::unetfc::VGNode> &pointer)
  const;

  /*! @brief Compute intersecting element and weight in tumor domain */
  void compute_elem_weights();

  /** @}*/

  /*! @brief Did we change network from last step */
  bool d_is_network_changed;

  /*! @brief Reference to model class */
  util::BaseModel *d_model_p;

  /*! @brief 1-D mesh */
  util::unetfc::ListStructure<util::unetfc::VGNode> VGM;

  /*! @brief System matrix for vessel pressure */
  gmm::row_matrix<gmm::wsvector<double>> A_VGM;

  /*! @brief System matrix for vessel nutrient */
  gmm::row_matrix<gmm::wsvector<double>> Ac_VGM;

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

  /*! @brief Old TAF concentration */
  std::vector<double> phi_TAF_old;

  double mu;

  double D_v;

  double D_v_3D;

  double D_TAF;

  double osmotic_sigma;

  std::string scenario;

  double h_3D;

  int N_3D;

  int N_tot_3D;

  unsigned int d_update_number;
};

} // namespace unetfc

} // namespace util

#endif // UTIL_UNETFC_NETWORK_H
