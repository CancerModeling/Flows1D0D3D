////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_MODEL_H
#define NETFV_MODEL_H

#include "umodel/model.hpp"
#include "modelUtil.hpp"
#include "systems/systems.hpp"
#include "unet/network.hpp"

// typedef network
typedef util::Network Net;

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace netfv {

void model_setup_run(int argc,
                     char **argv,
                     const std::string &filename,
                     Parallel::Communicator *comm);

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Model : public util::BaseModel {

public:
  /*! @brief Constructor */
  Model(int argc, char **argv,
        const std::string &filename,
        Parallel::Communicator *comm,
        InpDeck &input,
        ReplicatedMesh &mesh,
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

  const Net &get_network() const {
    return d_network;
  }
  Net &get_network() {
    return d_network;
  }

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

  /*! @brief Run model */
  void run() override ;

private:

  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step) override ;

  /*! @brief Solves tumor system */
  void solve_system() override ;

  /*! @brief Compute quantity of interest */
  void compute_qoi() override ;

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

  /*! @brief Network class */
  Net d_network;

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

  /*! @brief Ids of system to be used in logger */
  int d_nec_id;
  int d_tum_id;
  int d_nut_id;
  int d_hyp_id;
  int d_taf_id;
  int d_ecm_id;
  int d_mde_id;
  int d_pres_id;
  int d_grad_taf_id;
  int d_vel_id;
  int d_pres_1d_id;
  int d_nut_1d_id;
};

} // namespace netfv

#endif // NETFV_MODEL_H
