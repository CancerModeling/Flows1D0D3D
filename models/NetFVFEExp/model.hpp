////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFVFEEXP_MODEL_H
#define NETFVFEEXP_MODEL_H

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
namespace netfvfeexp {

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
      return d_pres;
    else if (system == "Nutrient")
      return d_nut;
    else if (system == "Prolific")
      return d_pro;
    else if (system == "Hypoxic")
      return d_hyp;
    else if (system == "Necrotic")
      return d_nec;
    else if (system == "TAF")
      return d_taf;
    else if (system == "TAF_Gradient")
      return d_grad_taf;
    else if (system == "ECM")
      return d_ecm;
    else if (system == "MDE")
      return d_mde;
    else if (system == "Velocity")
      return d_vel;
    else if (system == "Tumor")
      return d_tum;
    else
      libmesh_error_msg("Invalid system = " + system + " name");
  }
  PressureAssembly &get_pres_assembly() { return d_pres; }
  NutAssembly &get_nut_assembly() { return d_nut; }
  ProAssembly &get_pro_assembly() { return d_pro; }
  HypAssembly &get_hyp_assembly() { return d_hyp; }
  NecAssembly &get_nec_assembly() { return d_nec; }
  TafAssembly &get_taf_assembly() { return d_taf; }
  GradTafAssembly &get_grad_taf_assembly() { return d_grad_taf; }
  EcmAssembly &get_ecm_assembly() { return d_ecm; }
  MdeAssembly &get_mde_assembly() { return d_mde; }
  VelAssembly &get_vel_assembly() { return d_vel; }
  TumAssembly &get_tum_assembly() { return d_tum; }

  std::vector<util::BaseAssembly *> get_all_assembly() override {
    return {&d_tum, &d_nut, &d_pro, &d_hyp, &d_nec, &d_taf,
            &d_ecm, &d_mde, &d_pres, &d_grad_taf, &d_vel};
  }

  // list of systems to be solved inside nonlinear iteration loop
  std::vector<util::BaseAssembly *> get_nl_solve_assembly() {
    return {&d_nut, &d_pro, &d_hyp, &d_nec, &d_mde, &d_ecm};
  }

  /*! @brief Run model */
  void run() override;

private:
  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step) override;

  /*! @brief Solves tumor system */
  void solve_system() override;

  // all 3D systems explicit including 3D+1D nutrient
  // 1D-3D coupling in Nutrient and pressure is implicit in all cases
  void solve_system_explicit();

  // Nutrient is explicit ie out of the nonlineat iterations
  void solve_system_nutrient_explicit();

  // all 3D systems implicit
  void solve_system_implicit();


  /*! @brief Compute quantity of interest */
  void compute_qoi() override;

  /*! @brief Solves 1D-3D pressure system
   * At given time step, it solves for 1D and 3D pressure in a nonlinear loop
   * iteratively.
   */
  void solve_pressure();
  void solve_nutrient();

  /*! @brief Network class */
  Net d_network;

  /*! @brief Assembly objects */
  TumAssembly d_tum;
  NutAssembly d_nut;
  ProAssembly d_pro;
  HypAssembly d_hyp;
  NecAssembly d_nec;
  TafAssembly d_taf;
  EcmAssembly d_ecm;
  MdeAssembly d_mde;
  PressureAssembly d_pres;
  GradTafAssembly d_grad_taf;
  VelAssembly d_vel;

  /*! Ghosting functor so that PetSc does not give error when coupling dofs
   * of element with neighboring elements */
  util::GhostingFunctorFV d_ghosting_fv;

  UniquePtr<NumericVector<Number>> d_err_check_pro;
  UniquePtr<NumericVector<Number>> d_err_check_hyp;
  UniquePtr<NumericVector<Number>> d_err_check_nec;
  UniquePtr<NumericVector<Number>> d_err_check_nut;
  UniquePtr<NumericVector<Number>> d_err_check_mde;
  UniquePtr<NumericVector<Number>> d_err_check_ecm;
  UniquePtr<NumericVector<Number>> d_err_check_pres;
};

} // namespace netfvfeexp

#endif // NETFVFEEXP_MODEL_H
