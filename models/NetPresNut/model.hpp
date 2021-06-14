////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETPRESNUT_MODEL_H
#define NETPRESNUT_MODEL_H

#include "modelUtil.hpp"
#include "systems/systems.hpp"
#include "umodel/model.hpp"
#include "unet/network.hpp"
#include "unet/network_vtk_writer.h"
#include "usystem/ghosting_functor.hpp"

// typedef network
typedef util::unet::Network Net;

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace netpresnut {

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
        TransientLinearImplicitSystem &nut,
        TransientLinearImplicitSystem &pres,
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
    else
      libmesh_error_msg("Error: system = " + system + " does not exist");
  }
  PressureAssembly &get_pres_assembly() { return d_pres_assembly; }
  NutAssembly &get_nut_assembly() { return d_nut_assembly; }

  /*! @brief Run model */
  void run() override;

private:
  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step) override;

  /*! @brief Solves tumor system */
  void solve_system() override;

  /*! @brief Compute quantity of interest */
  void compute_qoi() override;

  /*! @brief Solves 1D-3D pressure system
   * At given time step, it solves for 1D and 3D pressure in a nonlinear loop
   * iteratively.
   */
  void solve_pressure();

  /*! @brief Network class */
  Net d_network;

  /*! @brief Saves the network as a vtk file. */
  util::unet::network_vtk_writer d_networkVtkWriter;
  util::unet::NetworkVTKWriterOld d_networkVtkWriterOld;

  /*! @brief Assembly objects */
  NutAssembly d_nut_assembly;
  PressureAssembly d_pres_assembly;

  /*! @brief Ids of system to be used in logger */
  int d_nut_id;
  int d_pres_id;
  int d_pres_1d_id;
  int d_nut_1d_id;

  /*! Ghosting functor so that PetSc does not give error when coupling dofs
   * of element with neighboring elements */
  util::GhostingFunctorFV d_ghosting_fv;
};

} // namespace netpresnut

#endif // NETPRESNUT_MODEL_H
