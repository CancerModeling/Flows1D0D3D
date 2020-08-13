////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TWOSP_MODEL_H
#define TWOSP_MODEL_H

#include "umodel/model.hpp"
#include "modelUtil.hpp"
#include "systems/systems.hpp"


/*!
 * @brief Namespace for simple two-species model
 */
namespace twosp {

void model_setup_run(int argc,
                     char **argv,
                     const std::string &filename,
                     Parallel::Communicator *comm);

/*!
 * @brief Two-species model consisting of tumor and nutrient
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
        TransientLinearImplicitSystem &tum,
        TransientLinearImplicitSystem &nut,
        util::Logger &log);

  /*! @brief Get various system classes */
  util::BaseAssembly &get_assembly(const std::string &system) override {
    if (system == "Nutrient")
      return d_nut_assembly;
    else if (system == "Tumor")
      return d_tum_assembly;
  }
  NutAssembly &get_nut_assembly() {return d_nut_assembly;}
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

  /*! @brief Assembly objects */
  TumAssembly d_tum_assembly;
  NutAssembly d_nut_assembly;

  /*! @brief Ids of system to be used in logger */
  int d_tum_id;
  int d_nut_id;
};

} // namespace twosp

#endif // TWOSP_MODEL_H
