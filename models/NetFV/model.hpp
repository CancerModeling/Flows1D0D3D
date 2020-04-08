////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_MODEL_H
#define NETFV_MODEL_H

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
namespace netfv {

void model_setup_run(int argc, char **argv, std::vector<double> &QOI_MASS,
                     const std::string &filename,
                     Parallel::Communicator *comm);

void create_mesh(netfv::InputDeck &input, ReplicatedMesh &mesh);

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Model {

public:
  /*! @brief Constructor */
  Model(int argc, char **argv, std::vector<double> &QOI_MASS,
        const std::string &filename, Parallel::Communicator *comm,
        netfv::InputDeck &input, ReplicatedMesh &mesh,
        EquationSystems &tum_sys,
        TransientLinearImplicitSystem &nec,
        TransientLinearImplicitSystem &tum,
        TransientLinearImplicitSystem &nut,
        TransientLinearImplicitSystem &hyp,
        TransientLinearImplicitSystem &taf,
        TransientLinearImplicitSystem &ecm,
        TransientLinearImplicitSystem &mde,
        TransientLinearImplicitSystem &pres,
        TransientLinearImplicitSystem &grad_taf,
        TransientLinearImplicitSystem &vel,
        util::Logger &log);

  /*! @brief Get equation system */
  const EquationSystems &get_system() const { return d_tum_sys; }
  EquationSystems &get_system() { return d_tum_sys; }

  const MeshBase &get_mesh() const { return d_mesh; }
  MeshBase &get_mesh() { return d_mesh; }

  const netfv::Network &get_network() const {
    return d_network;
  }
  netfv::Network &get_network() {
    return d_network;
  }

  const netfv::ListStructure<netfv::VGNode> &get_network_mesh() const {
    return d_network.get_mesh();
  }
  netfv::ListStructure<netfv::VGNode> &get_network_mesh() {
    return d_network.get_mesh();
  }

  /*! @brief Get input deck */
  const InputDeck &get_input_deck() const { return d_input; }
  InputDeck &get_input_deck() { return d_input; }

  /*! @brief Get network data */
  bool is_network_changed() const { return d_network.d_is_network_changed; }

  /*! @brief Get various system classes */
  PressureAssembly &get_pres_assembly() {return d_pres_assembly;}
  NutAssembly &get_nut_assembly() {return d_nut_assembly;}
  TumAssembly &get_tum_assembly() {return d_tum_assembly;}
  HypAssembly &get_hyp_assembly() {return d_hyp_assembly;}
  NecAssembly &get_nec_assembly() {return d_nec_assembly;}
  TafAssembly &get_taf_assembly() {return d_taf_assembly;}
  GradTafAssembly &get_grad_taf_assembly() {return d_grad_taf_assembly;}
  EcmAssembly &get_ecm_assembly() {return d_ecm_assembly;}
  MdeAssembly &get_mde_assembly() {return d_mde_assembly;}
  VelAssembly &get_vel_assembly() {return d_vel_assembly;}


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

  /*! @brief Current nonlinear step */
  unsigned int d_nonlinear_step;

  /*! @brief Is this output time step */
  bool d_is_output_step;

  /*! @brief Is this growth time step */
  bool d_is_growth_step;

  /*! @brief Reference to logger */
  util::Logger &d_log;

private:

  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step,
                    std::vector<double> *QOI_MASS = nullptr);

  /*! @brief Projects species concentrations to their physical range [0,1] */
  void project_solution_to_physical_range(const MeshBase &mesh,
                                          TransientLinearImplicitSystem &sys);

  /*! @brief Solves tumor system */
  void solve_system();

  /*! @brief Solves 1D-3D pressure system
   * At given time step, it solves for 1D and 3D pressure in a nonlinear loop
   * iteratively.
   */
  void solve_pressure();

  /*! @brief Solving sub-systems
   *
   * Description of subsystems
   *
   * ## test_nut
   *
   * - Nonlinear iterations to solve 1D and 3D nutrients
   * - Good to test the coupling of nutrient
   * - Maybe it is good to specify d_decouple_nutrients = false in input file
   *
   * - Good to test 1D/3D nutrient coupling
   *
   * ## test_nut_2
   *
   * - Solves 1D nutrient and then solves 3D nutrient
   * - This is recommended when d_decouple_nutrients = false
   *
   * - Good to test 1D/3D nutrient coupling
   *
   * ## test_taf
   *
   * - Solves only TAF equations
   * - Adds artificial cylindrical TAF source (along z axis)
   *    - specify center and radius of cylindrical source in input file
   *
   * - Good to test TAF equation and also growth of network
   *
   * ## test_taf_2
   *
   * - Solves pressure and TAF equations
   * - Adds artificial cylindrical TAF source (along z axis)
   *    - specify center and radius of cylindrical source in input file
   *
   * - Good to test TAF equation and also growth of network
   *
   * - Also test the pressure in growing network
   *
   * ## test_tum
   *
   * - Only solves tumor species
   * - Nutrient is constant
   * - Artificial source is added on cylinder along z-axis with
   *    - center = (L - 2R, L - 2R, 0)
   *    - Radius = 0.05 * L
   *    - L = size of domain
   *
   * - Good to test tumor species
   *
   * ## test_tum_2
   *
   * - Solves nutrient and tumor species
   * - Adds artificial cylindrical nutrient source (along z axis)
   *    - specify center and radius of cylindrical source in input file
   *
   * - Good to test coupling between nutrient and tumor species
   *
   * ## test_net_tum
   *
   * - Solves for 1D-3D pressure, 1D-3D nutrients, and all tumor species
   *
   * - Also solves for TAF at certain steps such as when we are producing
   * output or when we are growing the network
   *
   * - Good to test full model with/without network growth in which we do not
   * solve for MDE and ECM and avoid solving for TAF every time step
   *
   * ## test_net_tum_2
   *
   * - Same as test_net_tum, except that it solves for pressure at the
   * beginning and after that the pressure in tissue domain and network is fixed
   *
   * - Good to test the couling of 1D-3D nutrient coupling and coupling of
   * network and nutrient with tumor species
   *
   */
  void test_nut();
  void test_nut_2();
  void test_taf();
  void test_taf_2();
  void test_tum();
  void test_tum_2();
  void test_net_tum();
  void test_net_tum_2();

  /*! @brief To store input parameters */
  netfv::InputDeck &d_input;

  /*! @brief Pointer to communicator */
  Parallel::Communicator *d_comm_p;

  /*! @brief Store the network mesh */
  ReplicatedMesh &d_mesh;

  /*! @brief Store the 2d/3d tumor system */
  EquationSystems &d_tum_sys;

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
  VelAssembly d_vel_assembly;
};

} // namespace netfv

#endif // NETFV_MODEL_H
