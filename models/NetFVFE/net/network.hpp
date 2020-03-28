////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_NETWORK_H
#define NETFVFE_NETWORK_H

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
namespace netfvfe {

// forward declare Model class
class Model;

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Network {

public:
  /*! @brief Constructor */
  Network(netfvfe::Model *model)
      : d_is_network_changed(false), d_model_p(model) {}

  const netfvfe::ListStructure<netfvfe::VGNode> &get_mesh() const { return VGM; }

  netfvfe::ListStructure<netfvfe::VGNode> &get_mesh() { return VGM; }

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

  void writeDataToVTK_VGM();

  void writeDataToVTKTimeStep_VGM(int timeStep);

  /** @}*/

  /**
   * @name Solve
   */
  /**@{*/

  /*! @brief Solve 1-d system */
  void solve_system();
  void assembleVGMSystemForPressure();

  void solveVGMforPressure();

  void assembleVGMSystemForNutrient();

  void solveVGMforNutrient();

  /** @}*/

  /**
   * @name Update
   */
  /**@{*/

  void refine1DMesh();

  /*! @brief Update network */
  void update_network();

  /*! @brief Compute intersecting element and weight in tumor domain */
  void compute_elem_weights();

  /** @}*/

  /*! @brief Did we change network from last step */
  bool d_is_network_changed;

  /*! @brief Reference to model class */
  netfvfe::Model *d_model_p;

  /*! @brief 1-D mesh */
  netfvfe::ListStructure<netfvfe::VGNode> VGM;

  /*! @brief System matrix for vessel pressure */
  gmm::row_matrix<gmm::wsvector<double>> A_VGM;

  /*! @brief System matrix for vessel nutrient */
  gmm::row_matrix<gmm::wsvector<double>> Ac_VGM;

  /*! @brief System force for vessel pressure */
  std::vector<double> b;

  /*! @brief System force for vessel nutrient */
  std::vector<double> b_c;

  /*! @brief Current vessel pressure */
  std::vector<double> P_v;

  /*! @brief Current vessel nutrient */
  std::vector<double> C_v;
  std::vector<double> C_v_old;

  double mu;

  double D_v;

  double osmotic_sigma;

};

} // namespace netfvfe

#endif // NETFVFE_NETWORK_H
