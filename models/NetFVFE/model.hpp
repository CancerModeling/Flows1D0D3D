////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFE_MODEL_H
#define NETFVFE_MODEL_H

#include "umodel/model.hpp"
#include "usystem/ghosting_functor.hpp"
#include "modelUtil.hpp"
#include "systems/systems.hpp"
#include "unet/network.hpp"

// typedef network
typedef util::unet::Network Net;

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace netfvfe {

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
        TransientLinearImplicitSystem &pro,
        TransientLinearImplicitSystem &nut,
        TransientLinearImplicitSystem &hyp,
        TransientLinearImplicitSystem &taf,
        TransientLinearImplicitSystem &ecm,
        TransientLinearImplicitSystem &mde,
        TransientLinearImplicitSystem &pres,
        TransientLinearImplicitSystem &grad_taf,
        TransientLinearImplicitSystem &vel,
        TransientLinearImplicitSystem &tum,
        util::Logger &log);

  const Net &get_network() const {
    return d_network;
  }
  Net &get_network() {
    return d_network;
  }

  /*! @brief Get various system classes */
  util::BaseAssembly &get_assembly(const std::string &system) override {
    if (system == "Pressure")
      return d_pres_assembly;
    else if (system == "Nutrient")
      return d_nut_assembly;
    else if (system == "Prolific")
      return d_pro_assembly;
    else if (system == "Hypoxic")
      return d_hyp_assembly;
    else if (system == "Necrotic")
      return d_nec_assembly;
    else if (system == "TAF")
      return d_taf_assembly;
    else if (system == "TAF_Gradient")
      return d_grad_taf_assembly;
    else if (system == "ECM")
      return d_ecm_assembly;
    else if (system == "MDE")
      return d_mde_assembly;
    else if (system == "Velocity")
      return d_vel_assembly;
    else if (system == "Tumor")
      return d_tum_assembly;
    else
      libmesh_error_msg("Invalid system = " + system + " name");
  }
  PressureAssembly &get_pres_assembly() {return d_pres_assembly;}
  NutAssembly &get_nut_assembly() {return d_nut_assembly;}
  ProAssembly &get_pro_assembly() {return d_pro_assembly;}
  HypAssembly &get_hyp_assembly() {return d_hyp_assembly;}
  NecAssembly &get_nec_assembly() {return d_nec_assembly;}
  TafAssembly &get_taf_assembly() {return d_taf_assembly;}
  GradTafAssembly &get_grad_taf_assembly() {return d_grad_taf_assembly;}
  EcmAssembly &get_ecm_assembly() {return d_ecm_assembly;}
  MdeAssembly &get_mde_assembly() {return d_mde_assembly;}
  VelAssembly &get_vel_assembly() {return d_vel_assembly;}
  TumAssembly &get_tum_assembly() {return d_tum_assembly;}

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

  /*! @brief Network class */
  Net d_network;

  /*! @brief Assembly objects */
  NecAssembly d_nec_assembly;
  ProAssembly d_pro_assembly;
  NutAssembly d_nut_assembly;
  HypAssembly d_hyp_assembly;
  TafAssembly d_taf_assembly;
  EcmAssembly d_ecm_assembly;
  MdeAssembly d_mde_assembly;
  PressureAssembly d_pres_assembly;
  GradTafAssembly d_grad_taf_assembly;
  VelAssembly d_vel_assembly;
  TumAssembly d_tum_assembly;

  /*! @brief Ids of system to be used in logger */
  int d_tum_id;
  int d_nec_id;
  int d_pro_id;
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

  /*! Ghosting functor so that PetSc does not give error when coupling dofs
   * of element with neighboring elements */
  util::GhostingFunctorFV d_ghosting_fv;
};

} // namespace netfvfe

#endif // NETFVFE_MODEL_H
