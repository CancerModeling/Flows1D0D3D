////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNET_NETWORK_H
#define UTIL_UNET_NETWORK_H

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
namespace unet {

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Network {

public:
  /*! @brief Constructor */
  Network(util::BaseModel *model)
      : d_is_network_changed(false), d_model_p(model), d_update_number(0) {}

  const util::unet::ListStructure<util::unet::VGNode> &get_mesh() const { return
  VGM; }

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

  void writeDataToVTKTimeStep_VGM(int timeStep);

  /** @}*/

  /**
   * @name Solve
   */
  /**@{*/

  /*! @brief Solve 1-d system */
  void assembleVGMSystemForPressure(BaseAssembly &pres_sys);

  void solveVGMforPressure(BaseAssembly &pres_sys);

  void assembleVGMSystemForNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  void solveVGMforNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys);

  void update_old_concentration() { C_v_old = C_v; };

  std::vector<util::unet::ElemWeights>
  compute_elem_weights_at_node(std::shared_ptr<util::unet::VGNode> &pointer) const;

  /** @}*/

  /**
   * @name Update
   */
  /**@{*/

  void refine1DMesh();

  /*! @brief Update network */
  void update_network(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys);

  /*! @brief Compute intersecting element and weight in tumor domain */
  void compute_elem_weights();

  unsigned int markApicalGrowth(std::string growth_type, BaseAssembly &taf_sys);
  unsigned int processApicalGrowthTAF(BaseAssembly &taf_sys,
                                      BaseAssembly &grad_taf_sys);
  std::shared_ptr<util::unet::VGNode>
  check_new_node(const int &parent_index,
                 const std::vector<double> &parent_coord,
                 const std::vector<double> &child_coord, const double &dist_tol,
                 const double &domain_size, unsigned int &check_code);

  void add_new_node(std::shared_ptr<util::unet::VGNode> &pointer, const double &child_r,
                    const std::vector<double> &child_end_point);
  void add_new_node_at_existing_node(std::shared_ptr<util::unet::VGNode> &pointer,
                                     std::shared_ptr<util::unet::VGNode> &near_node);

  void get_taf_and_gradient(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys,
                            double &taf_val, Point &grad_taf_val,
                            const std::vector<double> &coord) const;

  /** @}*/

  /*! @brief Did we change network from last step */
  bool d_is_network_changed;

  /*! @brief Reference to model class */
  util::BaseModel *d_model_p;

  /*! @brief 1-D mesh */
  util::unet::ListStructure<util::unet::VGNode> VGM;

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

  unsigned int d_update_number;
};

} // namespace unet

} // namespace util

#endif // UTIL_UNET_NETWORK_H
