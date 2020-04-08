////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef AVAFV_MODEL_H
#define AVAFV_MODEL_H

#include "netinp/inp.hpp"
#include "systems/systems.hpp"
#include "utilLibs.hpp"
#include "utils.hpp"
#include <numeric>
#include <string>
#include <vector>

// typedef input deck so that change its namespace or name does not effect
// the rest of the code
typedef util::InputDeck InpDeck;

/*!
 * @brief Namespace for coupled 3d tumor growth model and 1d blood flow
 * network model. See
 * docs/NetTum/network_and_tumor_model.pdf for more details.
 */
namespace avafv {

void model_setup_run(int argc, char **argv, std::vector<double> &QOI_MASS,
                     const std::string &filename,
                     Parallel::Communicator *comm);

void create_mesh(InpDeck &input, ReplicatedMesh &mesh);

/*!
 * @brief Coupled 3D-1D tumor growth model driver
 */
class Model {

public:
  /*! @brief Constructor */
  Model(int argc, char **argv, std::vector<double> &QOI_MASS,
        const std::string &filename, Parallel::Communicator *comm,
        InpDeck &input, ReplicatedMesh &mesh,
        EquationSystems &tum_sys,
        TransientLinearImplicitSystem &nec,
        TransientLinearImplicitSystem &tum,
        TransientLinearImplicitSystem &nut,
        TransientLinearImplicitSystem &hyp,
        TransientLinearImplicitSystem &taf,
        TransientLinearImplicitSystem &grad_taf);

  /*! @brief Get equation system */
  const EquationSystems &get_system() const { return d_tum_sys; }
  EquationSystems &get_system() { return d_tum_sys; }

  const MeshBase &get_mesh() const { return d_mesh; }
  MeshBase &get_mesh() { return d_mesh; }

  /*! @brief Get input deck */
  const InpDeck &get_input_deck() const { return d_input; }
  InpDeck &get_input_deck() { return d_input; }

  /*! @brief Get various system classes */
  NutAssembly &get_nut_assembly() {return d_nut_assembly;}
  TumAssembly &get_tum_assembly() {return d_tum_assembly;}
  HypAssembly &get_hyp_assembly() {return d_hyp_assembly;}
  NecAssembly &get_nec_assembly() {return d_nec_assembly;}
  TafAssembly &get_taf_assembly() {return d_taf_assembly;}
  GradTafAssembly &get_grad_taf_assembly() {return d_grad_taf_assembly;}


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

private:

  /*! @brief Output results of tumor system and network system */
  void write_system(const unsigned int &t_step,
                    std::vector<double> *QOI_MASS = nullptr);

  /*! @brief Projects species concentrations to their physical range [0,1] */
  void project_solution_to_physical_range(const MeshBase &mesh,
                                          TransientLinearImplicitSystem &sys);

  /*! @brief Solves tumor system */
  void solve_system();

  /*! @brief Solves tumor system */
  void test_tum();
  void test_tum_2();

  /*! @brief To store input parameters */
  InpDeck &d_input;

  /*! @brief Pointer to communicator */
  Parallel::Communicator *d_comm_p;

  /*! @brief Store the network mesh */
  ReplicatedMesh &d_mesh;

  /*! @brief Store the 2d/3d tumor system */
  EquationSystems &d_tum_sys;

  /*! @brief Assembly objects */
  NecAssembly d_nec_assembly;
  TumAssembly d_tum_assembly;
  NutAssembly d_nut_assembly;
  HypAssembly d_hyp_assembly;
  TafAssembly d_taf_assembly;
  GradTafAssembly d_grad_taf_assembly;

  /*! @brief Test name */
  std::string d_test_name;
};

} // namespace avafv

#endif // AVAFV_MODEL_H
