////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_MODEL_H
#define NETFVFE_MODEL_H

#include "inp/inp.hpp"
#include "net/network.hpp"
#include "systems/systems.hpp"
#include "utilLibs.hpp"
#include "utils.hpp"
#include <numeric>
#include <string>
#include <vector>

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace netfvfe {

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Model {

public:
  /*! @brief Constructor */
  Model(int argc, char **argv, std::vector<double> &QOI_MASS,
        const std::string &filename, Parallel::Communicator *comm);

  /*! @brief Get equation system */
  const EquationSystems &get_system() const { return d_tum_sys; }
  EquationSystems &get_system() { return d_tum_sys; }

  const MeshBase &get_mesh() const { return d_mesh; }
  MeshBase &get_mesh() { return d_mesh; }

  const netfvfe::Network &get_network() const {
    return d_network;
  }
  netfvfe::Network &get_network() {
    return d_network;
  }

  const netfvfe::ListStructure<netfvfe::VGNode> &get_network_mesh() const {
    return d_network.get_mesh();
  }
  netfvfe::ListStructure<netfvfe::VGNode> &get_network_mesh() {
    return d_network.get_mesh();
  }

  /*! @brief Get input deck */
  const InputDeck &get_input_deck() const { return d_input; }
  InputDeck &get_input_deck() { return d_input; }

  /*! @brief Get network data */
  bool is_network_changed() const { return d_network.d_is_network_changed; }

public:
  /*! @brief To store input parameters */
  unsigned int d_step;

  /*! @brief Current time */
  double d_time;

  /*! @brief Current time step */
  double d_dt;

  /*! @brief hmax and hmin */
  double d_hmin;
  double d_hmax;

  /*! @brief Bounding box */
  std::pair<Point, Point> d_bounding_box;

private:
  /*! @brief Create mesh for 2d/3d pde and 1d network */
  void create_mesh();

  /*! @brief Setup 2d/3d tumor system and 1d network system */
  void setup_system();

  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step,
                    std::vector<double> *QOI_MASS = nullptr);

  /*! @brief Projects species concentrations to their physical range [0,1] */
  void project_solution_to_physical_range(const MeshBase &mesh,
                                          TransientLinearImplicitSystem &sys);

  /*! @brief Solves tumor system */
  void solve_system();

  /*! @brief Solves tumor system */
  void test_nutrient();
  void test_nutrient_2();
  void test_taf();
  void test_taf_2();

  /*! @brief Solves tumor system */
  void test_pressure();

  /*! @brief To store input parameters */
  netfvfe::InputDeck d_input;

  /*! @brief Pointer to communicator */
  Parallel::Communicator *d_comm_p;

  /*! @brief Store the network mesh */
  ReplicatedMesh d_mesh;

  /*! @brief Store the 2d/3d tumor system */
  EquationSystems d_tum_sys;

  /*! @brief Network class */
  Network d_network;

  /*! @brief Assembly objects */
  NecAssembly d_nec_assembly;
  TumAssembly d_tum_assembly;
  NutAssembly d_nut_assembly;
  HypAssembly d_hyp_assembly;
  TafAssembly d_taf_assembly;
  EcmAssembly d_ecm_assembly;
  MdeAssembly d_mde_assembly;
  PressureAssembly d_pres_assembly;
  GradTafAssembly d_grad_taf_assembly;
  VelocityAssembly d_vel_assembly;
};

} // namespace netfvfe

#endif // NETFVFE_MODEL_H
