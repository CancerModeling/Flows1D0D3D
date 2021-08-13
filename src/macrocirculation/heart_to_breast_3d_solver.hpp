////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_HEART_TO_BREAST_3_D_SOLVER_HPP
#define TUMORMODELS_HEART_TO_BREAST_3_D_SOLVER_HPP

#include "mpi.h"
#include "base_model.hpp"
#include "heart_to_breast_3d_systems.hpp"
#include "outlet_weight_functions.hpp"
#include <memory>

namespace macrocirculation {

/*! @brief Collects input parameters in this class for 3D system. */
struct HeartToBreast3DSolverInputDeck {
  HeartToBreast3DSolverInputDeck(const std::string &filename = "");

  /*! @brief Read parameter from input file. */
  void read_parameters(const std::string &filename);

  /*! @brief Print parameters for debug. */
  std::string print_str();

  /*! @brief Capillary hydraulic conductivity. */
  double d_K_cap;
  /*! @brief Tissue hydraulic conductivity. */
  double d_K_tis;
  /*! @brief Artery - Capillary exchange permeability. */
  double d_Lp_cap;
  /*! @brief Capillary - Tissue exchange permeability. */
  double d_Lp_tis;
  /*! @brief Final simulation time. */
  double d_T;
  /*! @brief Size of time step. */
  double d_dt;
  /*! @brief Mesh size in case we create uniform mesh. */
  double d_h;
  /*! @brief Name of mesh file in case we load mesh from file. */
  std::string d_mesh_file;
  /*! @brief Path to output results. */
  std::string d_out_dir;
};

class HeartToBreast3DSolver : public BaseModel {
public:
  HeartToBreast3DSolver(MPI_Comm mpi_comm,
                        lm::Parallel::Communicator *libmesh_comm,
                        HeartToBreast3DSolverInputDeck &input,
                        lm::ReplicatedMesh &mesh,
                        lm::EquationSystems &eq_sys,
                        lm::TransientLinearImplicitSystem &p_cap,
                        lm::TransientLinearImplicitSystem &p_tis,
                        lm::ExplicitSystem &K_cap_field,
                        lm::ExplicitSystem &Lp_cap_field,
                        Logger &log);

  // first define virtual functions
  void run() override {};
  void write_system(const unsigned int &t_step) override {};
  void solve_system() override {};
  void compute_qoi() override {};

  /*! @brief Setup artificial perfusion outlets . */
  void setup_random_outlets(unsigned int num_perf_outlets = 10);

  /*! @brief Perform secondary setup remaining after constructor. */
  void setup();

  /*! @brief Solve 3D system; currently this includes capillary and tissue pressures. */
  void solve() {};

  /*! @brief Output 3D system to a .pvtu file. */
  void write_output() {};

  /*! @brief Output perfusion data to a .vtu file. */
  void write_perfusion_output(std::string out_file);

  /*! @brief Return current time. */
  double get_time() const {return 0.; };

  /*! @brief Set output folder name. */
  void set_output_folder(std::string output_dir) {};

public:
  /*! @brief MPI comm. (Note that we have another comm from libmesh defined in BaseModel class) */
  MPI_Comm d_mpi_comm;

  /*! @brief Input parameter set. */
  HeartToBreast3DSolverInputDeck &d_input;

  /*! @brief Capillary pressure assembly. */
  CapillaryPressure d_p_cap;
  /*! @brief Tissue pressure assembly. */
  TissuePressure d_p_tis;
  /*! @brief Capillary hydraulic conductivity parameter (nonhomogeneous, may vary element-wise). */
  lm::ExplicitSystem &d_K_cap_field;
  /*! @brief Artery - Capillary exchange permeability parameter (nonhomogeneous, may vary element-wise). */
  lm::ExplicitSystem &d_Lp_cap_field;

  /*!
   * @brief Coordinates of perfusion outlets.
   * TODO: Change the perfusion outlet data from the 1D0D interface
   */
  std::vector<lm::Point> d_perf_pts;
  /*!
   * @brief Radius of arteries at the perfusion outlet.
   * TODO: Change the perfusion outlet data from the 1D0D interface
   */
  std::vector<double> d_perf_radii;
  /*!
   * @brief Radius of neighborhood for perfusion.
   * TODO: Change the perfusion outlet data from the 1D0D interface
   */
  std::vector<double> d_perf_ball_radii;
  /*!
   * @brief Perfusion outlet pressure vector.
   * TODO: Change the perfusion outlet data from the 1D0D interface
   */
  std::vector<double> d_perf_pres;
  /*! @brief Perfusion weight functions for each outlet. */
  std::vector<std::unique_ptr<BaseOutletRadial>> d_perf_fns;
  /*! @brief List of 3D mesh elements affected by perfusion from the outlet. */
  std::vector<std::vector<lm::dof_id_type>> d_perf_elems_3D;
  /*!
   * @brief Coefficient appearing in total perfused mass to the tissue from artery outlet.
   * NOTE: This should be returned to the 1D0D system
   */
  std::vector<double> d_perf_coeff_a;
  /*!
   * @brief Coefficient appearing in total perfused mass to the tissue from artery outlet.
   * NOTE: This should be returned to the 1D0D system
   */
  std::vector<double> d_perf_coeff_b;

};

} // namespace macrocirculation


#endif //TUMORMODELS_HEART_TO_BREAST_3_D_SOLVER_HPP