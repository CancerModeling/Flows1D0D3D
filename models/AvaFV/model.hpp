////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef AVAFV_MODEL_H
#define AVAFV_MODEL_H

#include "modelUtil.hpp"
#include "systems/systems.hpp"
#include "umodel/model.hpp"

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace avafv {

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
        TransientLinearImplicitSystem &grad_taf,
        util::Logger &log);

  /*! @brief Get various system classes */
  util::BaseAssembly &get_assembly(const std::string &system) {
    if (system == "Nutrient")
      return d_nut_assembly;
    else if (system == "Tumor")
      return d_tum_assembly;
    else if (system == "Hypoxic")
      return d_hyp_assembly;
    else if (system == "Necrotic")
      return d_nec_assembly;
    else if (system == "TAF")
      return d_taf_assembly;
    else if (system == "TAF_Gradient")
      return d_grad_taf_assembly;

    throw std::runtime_error( "unknown system " + system );
  }

  NutAssembly &get_nut_assembly() { return d_nut_assembly; }
  TumAssembly &get_tum_assembly() { return d_tum_assembly; }
  HypAssembly &get_hyp_assembly() { return d_hyp_assembly; }
  NecAssembly &get_nec_assembly() { return d_nec_assembly; }
  TafAssembly &get_taf_assembly() { return d_taf_assembly; }
  GradTafAssembly &get_grad_taf_assembly() { return d_grad_taf_assembly; }

  /*! @brief Run model */
  void run() override;

private:
  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step) override;

  /*! @brief Solves tumor system */
  void solve_system() override;

  /*! @brief Compute quantity of interest */
  void compute_qoi() override;

  /*! @brief Solving sub-systems
   *
   * Description of subsystems
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
   */
  void test_tum();
  void test_tum_2();

  /*! @brief Assembly objects */
  NecAssembly d_nec_assembly;
  TumAssembly d_tum_assembly;
  NutAssembly d_nut_assembly;
  HypAssembly d_hyp_assembly;
  TafAssembly d_taf_assembly;
  GradTafAssembly d_grad_taf_assembly;
};

} // namespace avafv

#endif // AVAFV_MODEL_H
