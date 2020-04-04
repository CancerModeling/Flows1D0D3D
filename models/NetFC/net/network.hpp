////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFC_NETWORK_H
#define NETFC_NETWORK_H

#include "utilLibs.hpp"
#include "../inp/inp.hpp"
#include "../systems/systems.hpp"
#include "gmm.h"
#include "utils.hpp"
#include "list_structure.hpp"
#include "nodes.hpp"
#include <string>
#include <vector>

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace netfc {

// forward declare Model class
class Model;

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Network {

public:
  /*! @brief Constructor */
  Network(netfc::Model *model)
      : d_is_network_changed(false), d_model_p(model), d_update_number(0) {}

  const netfc::ListStructure<netfc::VGNode> &get_mesh() const { return VGM; }

  netfc::ListStructure<netfc::VGNode> &get_mesh() { return VGM; }

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

  /** @}*/

  /**
   * @name Solve
   */
  /**@{*/

  /*! @brief Solve 3D1D-system */
  void assembleVGMSystemForNutrients();

  void assemble3D1DSystemForPressure();

  void assemble3D1DSystemForNutrients();

  void solveVGMforNutrient(  int timeStep, double time );

  void solve3D1DFlowProblem( int timeStep, double time );

  void solve3D1DNutrientProblem( int timeStep, double time );

  void solve3DProlificCellProblem( int timeStep, double time );

  /** @}*/

  /**
   * @name Update
   */
  /**@{*/

  void refine1DMesh();

  /*! @brief Compute intersecting element and weight in tumor domain */
  void compute_elem_weights();

  double getDirichletValue( std::vector< double > center_face, double L_p, double radius );

  double getK1D( double s, double L_p, double radius );

  void writeDataToVTK3D_Pressure(std::vector<double> P_3D, std::vector< std::vector<double> > V_3D, int N_3D, double h_3D, int timeStep);

  void writeDataToVTK3D_Nutrients(std::vector<double> C_3D, int N_3D, double h_3D, int timeStep);

  void rescaleSecombData( std::vector< std::vector<double> >& vertices, std::vector<double>& pressures, std::vector<double>& radii, double epsilon );

  void updateNetwork();

  void markApicalGrowth();

  unsigned int processApicalGrowth();

  /** @}*/

  /*! @brief Did we change network from last step */
  bool d_is_network_changed;

  /*! @brief Reference to model class */
  netfc::Model *d_model_p;

  /*! @brief 1-D mesh */
  netfc::ListStructure<netfc::VGNode> VGM;

  /*! @brief System matrix for vessel nutrient */
  gmm::row_matrix<gmm::wsvector<double>> Ac_VGM;

  /*! @brief System matrix for 3D1D pressure flow model*/
  gmm::row_matrix<gmm::wsvector<double>> A_flow_3D1D;

  /*! @brief System matrix for 3D1D nutrient model*/
  gmm::row_matrix<gmm::wsvector<double>> A_nut_3D1D;

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

  /*! @brief Current 3D1D nutrient concentration */
  std::vector<double> phi_sigma;

  /*! @brief Old 3D1D nutrient concentration */
  std::vector<double> phi_sigma_old;

  /*! @brief Current prolific cell concentration */
  std::vector<double> phi_P;

  /*! @brief Old prolific cell concentration */
  std::vector<double> phi_P_old;

  /*! @Current solution vector prolific cell concentration */
  std::vector<double> sol_Vec_Phi_P;

  /*! @Old solution vector prolific cell concentration */
  std::vector<double> sol_Vec_Phi_P_old;

  double mu;

  double D_v;

  double osmotic_sigma;

  unsigned int d_update_number;

  std::string scenario;

  double h_3D;

  int N_3D;

  int N_tot_3D;

};

} // namespace netfc

#endif // NETFC_NETWORK_H
