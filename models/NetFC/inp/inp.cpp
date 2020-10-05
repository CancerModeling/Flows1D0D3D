////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "inp.hpp"
#include "csv.hpp"
#include "libmesh/getpot.h"

void netfc::ModelDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_dim = input("dimension", 2);
  d_domain_params[0] = input("domain_xmin", 0.);
  d_domain_params[2] = input("domain_ymin", 0.);
  d_domain_params[1] = input("domain_xmax", 1.);
  d_domain_params[3] = input("domain_ymax", 1.);
  if (d_dim == 2)
    d_domain_params[4] = 0.;
  else if (d_dim == 3)
    d_domain_params[4] = input("domain_zmin", 0.);
  if (d_dim == 2)
    d_domain_params[5] = 0.;
  else if (d_dim == 3)
    d_domain_params[5] = input("domain_zmax", 1.);

  // perform check
  if (d_domain_params[1] < d_domain_params[0] + 1.0e-8) {
    libmesh_error_msg("Check domain size in input file");
    exit(1);
  }
  if (d_domain_params[3] < d_domain_params[2] + 1.0e-8) {
    libmesh_error_msg("Check domain size in input file");
    exit(1);
  }
  if (d_dim > 2) {
    if (d_domain_params[5] < d_domain_params[4] + 1.0e-8) {
      libmesh_error_msg("Check domain size in input file");
      exit(1);
    }
  }

  d_assembly_method = input("assembly_method", 1);
}

void netfc::ModelDeck::print(unsigned int level) {

  out << "#\n";
  out << "# ModelDeck\n";
  out << "#\n\n";

  out << "# Dimension of the domain \n";
  out << "# Default: 2\n";
  out << "dimension = " << d_dim << "\n\n";

  out << "# Domain starting x-coordinate \n";
  out << "# Default: 0\n";
  out << "domain_xmin = " << d_domain_params[0] << "\n\n";

  out << "# Domain ending x-coordinate \n";
  out << "# Default: 1\n";
  out << "domain_xmax = " << d_domain_params[1] << "\n\n";

  out << "# Domain starting y-coordinate \n";
  out << "# Default: 0\n";
  out << "domain_ymin = " << d_domain_params[2] << "\n\n";

  out << "# Domain ending y-coordinate \n";
  out << "# Default: 1\n";
  out << "domain_ymax = " << d_domain_params[3] << "\n\n";

  out << "# Domain starting z-coordinate \n";
  out << "# Default: 0\n";
  out << "domain_zmin = " << d_domain_params[4] << "\n\n";

  out << "# Domain ending z-coordinate \n";
  out << "# Default: 1\n";
  out << "domain_zmax = " << d_domain_params[5] << "\n\n";

  out << "# Assembly method to be used in creating Stiffness matrix and Force"
         " vector \n";
  out << "# Default: 2\n";
  out << "# Description: Depending on how the weak form is written down we "
         "can have many ways to assemble matrix and vector\n";
  out << "# 1- Assembly method as described in the draft\n";
  out << "# 2- Assembly method in which various species concentrations are "
         "projected to the physical range [0,1]\n";
  out << "# 3- Similar to 2 but now with every mass source term in the "
         "right-hand side of the weak form\n";
  out << "assembly_method = " << d_assembly_method << "\n\n";
}

void netfc::RestartDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_restart = input("restart", false);
  d_mesh_restart_file = input("mesh_restart_file", "");
  d_sol_restart_file = input("solution_restart_file", "");
}

void netfc::RestartDeck::print(unsigned int level) {

  out << "#\n";
  out << "# RestartDeck\n";
  out << "#\n\n";

  out << "# Restart flag - Are we restarting simulation \n";
  out << "# Default: false\n";
  out << "restart = " << d_restart << "\n\n";

  out << "# Mesh file for restart \n";
  out << "# Default: <empty string>\n";
  out << "mesh_restart_file = " << d_mesh_restart_file << "\n\n";

  out << "# Solution file for restart \n";
  out << "# Default: <empty string>\n";
  out << "solution_restart_file = " << d_sol_restart_file << "\n\n";
}

void netfc::MeshDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_mesh_filename = input("mesh_name", "");
  d_num_elems = input("mesh_n_elements", 0);
  if (d_num_elems == 0) {
    d_num_elems_vec[0] = input("mesh_n_elements_x", 50);
    d_num_elems_vec[1] = input("mesh_n_elements_y", 50);
    d_num_elems_vec[2] = input("mesh_n_elements_z", 50);
  } else {
    d_num_elems_vec[0] = d_num_elems;
    d_num_elems_vec[1] = d_num_elems;
    d_num_elems_vec[2] = d_num_elems;
  }
  d_n_global_refinement = input("mesh_n_global_refinement", 0);
  d_read_mesh_flag = input("mesh_read_file", false);
}

void netfc::MeshDeck::print(unsigned int level) {

  out << "#\n";
  out << "# MeshDeck\n";
  out << "#\n\n";

  out << "# Mesh filename (if reading mesh from file) \n";
  out << "# Default: <empty string>\n";
  out << "mesh_name = " << d_mesh_filename << "\n\n";

  out << "# Number of elements for discretization in each coordinate \n";
  out << "# Default: 100\n";
  out << "mesh_n_elements = " << d_num_elems << "\n\n";

  out << "# Global refinement level \n";
  out << "# Default: 0\n";
  out << "mesh_n_global_refinement = " << d_n_global_refinement << "\n\n";

  out << "# Read mesh from file? If true, must specify mesh file using tag "
         "mesh_name above\n";
  out << "# Default: false\n";
  out << "mesh_read_file = " << d_read_mesh_flag << "\n\n";
}

void netfc::TimeDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_dt = input("time_step", 0.05);
  d_init_time = input("initial_time", 0.);
  d_final_time = input("final_time", 1.);
  d_max_time_steps = input("max_time_steps", 10000);

  d_steps = (d_final_time - d_init_time) / d_dt;

  d_init_step = input("initial_step", 0);
}

void netfc::TimeDeck::print(unsigned int level) {

  out << "#\n";
  out << "# TimeDeck\n";
  out << "#\n\n";

  out << "# Size of time step \n";
  out << "# Default: 0.05\n";
  out << "time_step = " << d_dt << "\n\n";

  out << "# Initial time of simulation \n";
  out << "# Default: 0\n";
  out << "# Description: In case of restart, this should be set to the time "
         "from which we are restarting the simulation. \n";
  out << "initial_time = " << d_init_time << "\n\n";

  out << "# Final time of simulation \n";
  out << "# Default: 1\n";
  out << "final_time = " << d_final_time << "\n\n";

  out << "# Maximum number of time steps \n";
  out << "# Default: 10000\n";
  out << "max_time_steps = " << d_max_time_steps << "\n\n";

  out << "# Initial step (only if restarting the simulation) \n";
  out << "# Default: 0\n";
  out << "# Description: When restarting the simulation, specify the time "
         "step from which we are restarting. This helps in creating "
         "output with compatible time tags\n";
  out << "initial_step = " << d_init_step << "\n\n";
}

void netfc::OutputDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_perform_output = input("perform_output", true);
  d_dt_output_interval = input("output_interval", 1);

  d_restart_save = input("restart_save", false);
  d_dt_restart_save_interval = input("restart_save_interval", 1);
}

void netfc::OutputDeck::print(unsigned int level) {

  out << "#\n";
  out << "# OutputDeck\n";
  out << "#\n\n";

  out << "# Perform output of simulation results or not?\n";
  out << "# Default: true\n";
  out << "perform_output = " << d_perform_output << "\n\n";

  out << "# Output interval for writing solution data to file \n";
  out << "# Default: 1\n";
  out << "# Description: Every N time steps the code writes simulation data "
         "to file.\n";
  out << "output_interval = " << d_dt_output_interval << "\n\n";

  out << "# Restart save flag \n";
  out << "# Default: false\n";
  out << "# Description: Set this flag as true if want to save the "
         "simulation data at every few time steps to a unique file. This "
         "saved data can then be used to restart the simulation.\n";
  out << "restart_save = " << d_restart_save << "\n\n";

  out << "# Restart save interval (effective only when restart_save is set to"
         " true) \n";
  out << "# Default: 1\n";
  out << "# Description: Specify time step interval for saving files for "
         "restart.\n";
  out << "restart_save_interval = " << d_dt_restart_save_interval << "\n\n";
}

void netfc::SolverDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_linear_max_iters = input("linear_solver_max_iter", 200);
  d_linear_tol = input("linear_solver_tol", 1.e-6);

  d_nonlin_max_iters = input("nonlinear_solver_max_iter", 200);
  d_nonlin_tol = input("nonlinear_solver_tol", 1.e-6);

  d_project_solution_to_physical_range =
    input("project_solution_to_phyiscal_range", false);
}

void netfc::SolverDeck::print(unsigned int level) {

  out << "#\n";
  out << "# SolverDeck\n";
  out << "#\n\n";

  out << "# Maximum number of iterations for linear solver\n";
  out << "# Default: 200\n";
  out << "linear_solver_max_iter = " << d_linear_max_iters << "\n\n";

  out << "# Tolerance for linear solver \n";
  out << "# Default: 1.e-6\n";
  out << "linear_solver_tol = " << d_linear_tol << "\n\n";

  out << "# Maximum number of iterations for nonlinear solver\n";
  out << "# Default: 200\n";
  out << "nonlinear_solver_max_iter = " << d_nonlin_max_iters << "\n\n";

  out << "# Tolerance for nonlinear solver \n";
  out << "# Default: 1.e-6\n";
  out << "nonlinear_solver_tol = " << d_nonlin_tol << "\n\n";

  out << "# Should we project solutions to physical range? (For now this does"
         " not seem to work so set it to false) \n";
  out << "# Default: false\n";
  out << "project_solution_to_phyiscal_range = "
      << d_project_solution_to_physical_range << "\n\n";
}

void netfc::NutrientDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_lambda_P = input("lambda_P", 0.5);
  d_lambda_A = input("lambda_A", 0.1);
  d_lambda_Ph = input("lambda_Ph", 0.1);
  d_D_sigma = input("D_sigma", 1.);
  d_delta_sigma = input("delta_sigma", 0.01);
  d_chi_c = input("chi_c", 0.035);
}

void netfc::NutrientDeck::print(unsigned int level) {

  out << "#\n";
  out << "# NutrientDeck\n";
  out << "#\n\n";

  out << "# Rate of nutrient consumption by prolific tumor cells\n";
  out << "# Default: 0.5\n";
  out << "Description: Negative source of nutrient mass\n";
  out << "lambda_P = " << d_lambda_P << "\n\n";

  out << "# Rate of apotosis of tumor cells\n";
  out << "# Default: 0.1\n";
  out << "Description: Positive source of nutrient mass as apotosis of cells "
         "results in increase of nutrients\n";
  out << "lambda_A = " << d_lambda_A << "\n\n";

  out << "# Rate of nutrient consumption by hypoxic tumor cells\n";
  out << "# Default: 0.1\n";
  out << "Description: Negative source of nutrient mass\n";
  out << "lambda_Ph = " << d_lambda_Ph << "\n\n";

  out << "# Diffusivity constant for nutrient species\n";
  out << "# Default: 1\n";
  out << "D_sigma = " << d_D_sigma << "\n\n";

  out << "# delta_sigma Constant\n";
  out << "# Default: 0.01\n";
  out << "delta_sigma = " << d_delta_sigma << "\n\n";

  out << "# Chemotaxis constant\n";
  out << "# Default: 0.035\n";
  out << "chi_c = " << d_chi_c << "\n\n";
}

void netfc::TumorDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_bar_M_P = input("bar_M_P", 200.);
  d_bar_E_phi_T = input("bar_E_phi_T", 0.045);
  d_epsilon_T = input("epsilon_T", 0.005);
}

void netfc::TumorDeck::print(unsigned int level) {

  out << "#\n";
  out << "# TumorDeck\n";
  out << "#\n\n";

  out << "# Constant for mobility of prolific tumor cells\n";
  out << "# Default: 200\n";
  out << "bar_M_P = " << d_bar_M_P << "\n\n";

  out << "# Constant associated to double-well potential of tumor cells \n";
  out << "# Default: 0.045\n";
  out << "bar_E_phi_T = " << d_bar_E_phi_T << "\n\n";

  out << "# Constant associated to interfacial energy of tumor cells \n";
  out << "# Default: 0.005\n";
  out << "epsilon_T = " << d_epsilon_T << "\n\n";
}

void netfc::HypoxicDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_bar_M_H = input("bar_M_H", 20.);
  d_lambda_HP = input("lambda_HP", 0.5);
  d_lambda_PH = input("lambda_PH", 0.5);
  d_lambda_HN = input("lambda_HN", 0.5);
  d_sigma_PH = input("sigma_PH", 0.4);
  d_sigma_HP = input("sigma_HP", 0.5);
  d_sigma_HN = input("sigma_HN", 0.2);
}

void netfc::HypoxicDeck::print(unsigned int level) {

  out << "#\n";
  out << "# HypoxicDeck\n";
  out << "#\n\n";

  out << "# Constant for mobility of hypoxic tumor cells\n";
  out << "# Default: 20\n";
  out << "bar_M_H = " << d_bar_M_H << "\n\n";

  out << "# Rate of conversion from hypoxic to prolific tumor cells\n";
  out << "# Default: 0.5\n";
  out << "Description: Negative source\n";
  out << "lambda_HP = " << d_lambda_HP << "\n\n";

  out << "# Rate of conversion from prolific to hypoxic tumor cells\n";
  out << "# Default: 0.5\n";
  out << "Description: Positive source\n";
  out << "lambda_PH = " << d_lambda_PH << "\n\n";

  out << "# Rate of conversion from hypoxic to necrotic tumor cells\n";
  out << "# Default: 0.5\n";
  out << "Description: Negative source\n";
  out << "lambda_HN = " << d_lambda_HN << "\n\n";

  out << "# Threshold of nutrient to trigger conversion from hypoxic to "
         "prolific tumor cells\n";
  out << "# Default: 0.4\n";
  out << "sigma_PH = " << d_sigma_PH << "\n\n";

  out << "# Threshold of nutrient to trigger conversion from prolific "
         "to hypoxic tumor cells\n";
  out << "# Default: 0.5\n";
  out << "sigma_HP = " << d_sigma_HP << "\n\n";

  out << "# Threshold of nutrient to trigger conversion from hypoxic to "
         "necrotic tumor cells\n";
  out << "# Default: 0.2\n";
  out << "sigma_HN = " << d_sigma_HN << "\n\n";
}

void netfc::NecroticDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_bar_M_N = input("bar_M_N", 0.);
}

void netfc::NecroticDeck::print(unsigned int level) {

  out << "#\n";
  out << "# NecroticDeck\n";
  out << "#\n\n";

  out << "# Constant for mobility of necrotic tumor cells\n";
  out << "# Default: 0\n";
  out << "bar_M_N = " << d_bar_M_N << "\n\n";
}

void netfc::TAFDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_D_TAF = input("D_TAF", 10.);
  d_delta_TAF = input("delta_TAF", 1.);
  d_lambda_TAF = input("lambda_TAF", 10.);
}

void netfc::TAFDeck::print(unsigned int level) {

  out << "#\n";
  out << "# TAFDeck\n";
  out << "#\n\n";

  out << "# Diffusivity constant for TAF species\n";
  out << "# Default: 10\n";
  out << "D_TAF = " << d_D_TAF << "\n\n";

  out << "# delta_TAF Constant\n";
  out << "# Default: 1\n";
  out << "delta_TAF = " << d_delta_TAF << "\n\n";

  out << "# Rate of TAF production by hpoxic cells\n";
  out << "# Default: 10\n";
  out << "lambda_TAF = " << d_lambda_TAF << "\n\n";
}

void netfc::ECMDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_lambda_ECM_D = input("lambda_ECM_D", 1.);
  d_lambda_ECM_P = input("lambda_ECM_P", 1.);
  d_bar_phi_ECM_P = input("bar_phi_ECM_P", 0.5);
  d_chi_h = input("chi_h", 0.035);

  // read ic information for ecm

  // locate the domain of stroma
  d_ecm_ic_data.d_type = input("ECM_ic_domain_type", "");

  // get the initial value of ecm in stroma
  d_ecm_ic_data.d_val = input("ECM_ic_val", 1.);

  // get the geometrical parameters
  unsigned int num_params = input("ECM_ic_num_params", 0);

  for (unsigned int i = 0; i < num_params; i++)
    d_ecm_ic_data.d_geom_params.push_back(input("ECM_ic_params_" + std::to_string(i + 1), 0.));
}

void netfc::ECMDeck::print(unsigned int level) {

  out << "#\n";
  out << "# ECMDeck\n";
  out << "#\n\n";

  out << "# Decay rate of ECM (fibronectin) due to MDE\n";
  out << "# Default: 1\n";
  out << "lambda_ECM_D = " << d_lambda_ECM_D << "\n\n";

  out << "# Production rate of ECM due to nutrients\n";
  out << "# Default: 1\n";
  out << "d_lambda_ECM_P = " << d_lambda_ECM_P << "\n\n";

  out << "# Threshold on ECM concentration for activation of production of "
         "ECM\n";
  out << "# Default: 0.5\n";
  out << "bar_phi_ECM_P = " << d_bar_phi_ECM_P << "\n\n";

  out << "# Haptotaxis factor for endothelial cells\n";
  out << "# Default: 0.035\n";
  out << "chi_h = " << d_chi_h << "\n\n";

  out << "# Type of domain where ECM is present initially (Stroma) \n";
  out << "# Default: none\n";
  out << "ECM_ic_domain_type = " << d_ecm_ic_data.d_type << "\n\n";

  out << "# Initial value of ECM in the domain \n";
  out << "# Default: 1\n";
  out << "ECM_ic_val = " << d_ecm_ic_data.d_val << "\n\n";

  out << "# Parameters which describe the geometry of ECM domain \n";
  out << "# Default: none\n";
  out << "ECM_ic_num_params = " << d_ecm_ic_data.d_geom_params.size() << "\n";
  for (unsigned int i = 0; i < d_ecm_ic_data.d_geom_params.size(); i++)
    out << "ECM_ic_params_" + std::to_string(i + 1)
        << d_ecm_ic_data.d_geom_params[i] << "\n";
  out << "\n";
}

void netfc::MDEDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_D_MDE = input("D_MDE", 10.);
  d_delta_MDE = input("delta_MDE", 1.);
  d_lambda_MDE_D = input("lambda_MDE_D", 1.);
  d_lambda_MDE_P = input("lambda_MDE_P", 1.);

  d_mde_ic_val = input("MDE_ic_val", 0.5);
}

void netfc::MDEDeck::print(unsigned int level) {

  out << "#\n";
  out << "# MDEDeck\n";
  out << "#\n\n";

  out << "# Diffusivity constant for MDE species\n";
  out << "# Default: 10\n";
  out << "D_MDE = " << d_D_MDE << "\n\n";

  out << "# delta_MDE Constant\n";
  out << "# Default: 1\n";
  out << "delta_MDE = " << d_delta_MDE << "\n\n";

  out << "# Natural decay rate of MDE \n";
  out << "# Default: 1\n";
  out << "lambda_MDE_D = " << d_lambda_MDE_D << "\n\n";

  out << "# Production rate of MDE due to hypoxic cells\n";
  out << "# Default: 1\n";
  out << "lambda_MDE_P = " << d_lambda_MDE_P << "\n\n";

  out << "# Initial value of MDE (only in region where initial condition on "
         "hypoxic cells are prescribed) \n";
  out << "# Default: 0.5\n";
  out << "MDE_ic_val = " << d_mde_ic_val << "\n\n";
}

void netfc::NutrientICDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_nut_ic_value = input("ic_nutrient_value", 0.5);
  // out << "nutrient ic = " << d_nut_ic_value;
}

void netfc::NutrientICDeck::print(unsigned int level) {

  out << "#\n";
  out << "# NutrientICDeck\n";
  out << "#\n\n";

  out << "# Initial value of nutrients \n";
  out << "# Default: 0.5\n";
  out << "ic_nutrient_value = " << d_nut_ic_value << "\n\n";
}

void netfc::TumorICDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  int dim = input("dimension", 2);

  // see if we read the ic data from csv file or from the input file
  bool read_from_csv = false;
  std::string csv_file = input("ic_tumor_file", "");
  if (!csv_file.empty())
    read_from_csv = true;

  if (read_from_csv) {

    // header
    // type (int), center (x,y,z), tum_radius(r1, r2, r3), hyp_radius(r1, r2,
    // r3)
    io::CSVReader<10> in(csv_file);

    in.read_header(io::ignore_extra_column, "type", "cx", "cy", "cz",
                   "tum_rx", "tum_ry", "tum_rz", "hyp_rx", "hyp_ry", "hyp_rz");

    int type;
    double cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz;
    while (in.read_row(type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx,
                       hyp_ry, hyp_rz)) {

      std::string ic_type;
      if (type == 1)
        ic_type = "tumor_spherical";
      else if (type == 2)
        ic_type = "tumor_elliptical";
      else if (type == 3)
        ic_type = "tumor_hypoxic_spherical";
      else if (type == 4)
        ic_type = "tumor_hypoxic_elliptical";

      auto data = TumorICData(ic_type, {cx, cy, cz}, {tum_rx, tum_ry, tum_rz},
                              {hyp_rx, hyp_ry, hyp_rz});

      d_tum_ic_data.push_back(data);
    }
  } else {

    // get how many tumor core are there
    int num_ic = input("ic_tumor_number", 0);
    d_tum_ic_data.resize(num_ic);

    for (unsigned int i = 0; i < num_ic; i++) {

      std::string postfix = "_" + std::to_string(i);
      if (num_ic == 1)
        postfix = "";

      auto data = TumorICData();

      data.d_ic_type = input("ic_tumor_type" + postfix, "");
      if (data.d_ic_type == "tumor_spherical" or
          data.d_ic_type == "tumor_hypoxic_spherical") {

        data.d_tum_ic_radius[0] = input("ic_tumor_radius" + postfix, 0.);

        if (data.d_ic_type == "tumor_hypoxic_spherical")
          data.d_hyp_ic_radius[0] = input("ic_hypoxic_radius" + postfix, 0.);

        if (data.d_hyp_ic_radius[0] < data.d_tum_ic_radius[0]) {
          libmesh_error_msg(
            "Error: Radius for hypoxic ic can not be smaller than"
            " tumor ic.");
          exit(1);
        }
      } else if (data.d_ic_type == "tumor_elliptical" or
                 data.d_ic_type == "tumor_hypoxic_elliptical") {

        data.d_tum_ic_radius[0] = input("ic_tumor_radius_x" + postfix, 0.);
        data.d_tum_ic_radius[1] = input("ic_tumor_radius_y" + postfix, 0.);
        data.d_tum_ic_radius[2] = input("ic_tumor_radius_z" + postfix, 0.);

        if (data.d_ic_type == "tumor_hypoxic_elliptical") {

          data.d_hyp_ic_radius[0] = input("ic_hypoxic_radius_x" + postfix, 0.);
          data.d_hyp_ic_radius[1] = input("ic_hypoxic_radius_y" + postfix, 0.);
          data.d_hyp_ic_radius[2] = input("ic_hypoxic_radius_z" + postfix, 0.);
        }

        if (data.d_hyp_ic_radius[0] < data.d_tum_ic_radius[0] or
            data.d_hyp_ic_radius[1] < data.d_tum_ic_radius[1] ||
            data.d_hyp_ic_radius[2] < data.d_tum_ic_radius[2]) {
          libmesh_error_msg(
            "Error: Radius for hypoxic ic can not be smaller than"
            " tumor ic.");
          exit(1);
        }
      }
      data.d_ic_center[0] = input("ic_tumor_center_x" + postfix, 0.);
      data.d_ic_center[1] = input("ic_tumor_center_y" + postfix, 0.);
      data.d_ic_center[2] = input("ic_tumor_center_z" + postfix, 0.);
      if (dim == 2)
        data.d_ic_center[2] = 0.;

      // add data
      d_tum_ic_data[i] = data;
    }
  }
}

void netfc::TumorICDeck::print(unsigned int level) {

  out << "#\n";
  out << "# TumorICDeck\n";
  out << "#\n\n";

  out << "# Tumor core type \n";
  out << "# Default: 0\n";
  out << "# Description: Following types of core has been implemented\n";
  out << "# 0- Spherical/circular core\n";
  //  out << "ic_tumor_type = " << d_tum_ic_type << "\n\n";

  out << "# Tumor core radius (if ic_tumor_type is 0)  \n";
  out << "# Default: 1\n";
  //  out << "ic_tumor_radius = " << d_tum_ic_radius[0] << "\n\n";

  out << "# Tumor core x-radius (if ic_tumor_type is 1)  \n";
  out << "# Default: 1\n";
  out << "# Description: If core is ellipsoidal then this is size of axis in "
         "x-direction \n";
  //  out << "ic_tumor_radius_x = " << d_tum_ic_radius[0] << "\n\n";

  out << "# Tumor core y-radius (if ic_tumor_type is 1)  \n";
  out << "# Default: 1\n";
  out << "# Description: If core is ellipsoidal then this is size of axis in "
         "y-direction \n";
  //  out << "ic_tumor_radius_y = " << d_tum_ic_radius[1] << "\n\n";

  out << "# Tumor core z-radius (if ic_tumor_type is 1)  \n";
  out << "# Default: 1\n";
  out << "# Description: If core is ellipsoidal then this is size of axis in "
         "z-direction \n";
  //  out << "ic_tumor_radius_z = " << d_tum_ic_radius[2] << "\n\n";

  out << "# x-coordinate of center of tumor core \n";
  out << "# Default: 0\n";
  //  out << "ic_tumor_center_x = " << d_tum_ic_center[0] << "\n\n";

  out << "# y-coordinate of center of tumor core \n";
  out << "# Default: 0\n";
  //  out << "ic_tumor_center_y = " << d_tum_ic_center[1] << "\n\n";

  out << "# z-coordinate of center of tumor core \n";
  out << "# Default: 0\n";
  //  out << "ic_tumor_center_z = " << d_tum_ic_center[2] << "\n\n";
}

void netfc::NutrientBCDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_nutrient_bc_north = input("bc_nutrient_north", false);
  d_nutrient_bc_south = input("bc_nutrient_south", false);
  d_nutrient_bc_east = input("bc_nutrient_east", false);
  d_nutrient_bc_west = input("bc_nutrient_west", false);
}

void netfc::NutrientBCDeck::print(unsigned int level) {

  out << "#\n";
  out << "# NutrientBCDeck\n";
  out << "#\n\n";

  out << "# Nutrient bc in north? \n";
  out << "# Default: false\n";
  out << "bc_nutrient_north = " << d_nutrient_bc_north << "\n\n";

  out << "# Nutrient bc in south? \n";
  out << "# Default: false\n";
  out << "bc_nutrient_south = " << d_nutrient_bc_south << "\n\n";

  out << "# Nutrient bc in east? \n";
  out << "# Default: false\n";
  out << "bc_nutrient_east = " << d_nutrient_bc_east << "\n\n";

  out << "# Nutrient bc in west? \n";
  out << "# Default: false\n";
  out << "bc_nutrient_west = " << d_nutrient_bc_west << "\n\n";
}

void netfc::NetworkDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  network_active = input("is_network_active", false);

  d_network_init_file = input("network_init_file", "");
  d_network_init_refinement = input("network_init_refinement", 1);
  d_num_points_length = input("network_discret_cyl_length", 2);
  d_num_points_angle = input("network_discret_cyl_angle", 2);
  d_coupling_method_theta = input("network_coupling_method_theta", 0.5);

  d_coupling_factor_p_t = input("coupling_factor_p_t", 1.);

  d_identify_vein_pres = input("identify_vein_pressure", 0.);

  d_net_direction_lambda_g = input("vessel_lambda_g", 0.);
  d_net_length_R_factor = input("vessel_R_factor", 0.);

  d_network_update_interval = input("network_update_interval", 1);

  d_log_normal_mean = input("log_normal_mean", 0.);
  d_log_normal_std_dev = input("log_normal_std_dev", 0.);

  d_net_radius_exponent_gamma = input("network_radius_exponent_gamma", 1.);

  d_no_branch_dist = input("network_no_branch_dist", 1);

  d_new_vessel_max_angle = input("network_new_veesel_max_angle", M_PI / 4.);

  d_branch_angle = input("network_branch_angle", M_PI / 8.);

  d_network_update_taf_threshold = input("network_update_taf_threshold", 0.);

  d_vessel_no_taf_effect_dist = input("network_vessel_no_taf_dist", 5);

  d_nonlocal_direction_search_num_points = input("network_nonlocal_search_num_points", 2);
  d_nonlocal_direction_search_length =
    input("network_nonlocal_search_length_factor", 10.);
  d_network_local_search = input("network_local_search", true);

  d_network_no_new_node_search_factor = input("network_no_new_node_search_factor", 0.5);
}

void netfc::NetworkDeck::print(unsigned int level) {
}

void netfc::Flow1DDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_init_vessel_mu = input("init_vessel_viscosity", 1.);
  d_in_pressure = input("vessel_in_pressure", 0.);
  d_in_nutrient = input("vessel_in_nutrient", 0.);
  d_in_nutrient_vein = input("vessel_in_nutrient_vein", 0.);
  d_blood_density = input("vessel_blood_density", 1.);
  d_D_sigma_v = input("vessel_D_sigma", 1.);

  d_osmotic_sigma = input("osmotic_reflection_coeff", 0.1);

  d_scenario = input("scenario", "scenario not available");
}

void netfc::Flow1DDeck::print(unsigned int level) {

  out << "#\n";
  out << "# Flow1DDeck\n";
  out << "#\n\n";

  out << "# Initial blood viscosity \n";
  out << "# Default: 1\n";
  out << "init_vessel_viscosity = " << d_init_vessel_mu << "\n\n";

  out << "# Boundary condition for pressure in network (In-pressure) \n";
  out << "# Default: 0\n";
  out << "vessel_in_pressure = " << d_in_pressure << "\n\n";
}

void netfc::FlowDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_tissue_flow_mu = input("tissue_flow_viscosity", 1.);
  d_tissue_flow_K = input("tissue_flow_K", 1.);
  d_tissue_flow_rho = input("tissue_flow_density", 1.);
  d_tissue_flow_L_p = input("tissue_flow_L_p", 1.);
  d_tissue_nut_L_s = input("tissue_nut_L_s", 1.);

  d_pressure_bc_val = input("tissue_pressure_bc_val", 0.);
  d_pressure_ic_val = input("tissue_pressure_ic_val", 0.);

  d_pressure_bc_north = input("bc_tissue_pressure_north", false);
  d_pressure_bc_south = input("bc_tissue_pressure_south", false);
  d_pressure_bc_east = input("bc_tissue_pressure_east", false);
  d_pressure_bc_west = input("bc_tissue_pressure_west", false);

  d_N_newton = input("N_newton", 0);
  d_omega = input("omega", 0.0);
}

void netfc::FlowDeck::print(unsigned int level) {

  out << "#\n";
  out << "# FlowDeck\n";
  out << "#\n\n";

  out << "# Viscosity of interstitial fluid in tissue \n";
  out << "# Default: 1\n";
  out << "tissue_flow_viscosity = " << d_tissue_flow_mu << "\n\n";

  out << "# Permeability constant of interstitial fluid in tissue \n";
  out << "# Default: 1\n";
  out << "tissue_flow_K = " << d_tissue_flow_K << "\n\n";

  out << "# Tissue pressure dirichlet boundary condition value \n";
  out << "# Default: 0\n";
  out << "tissue_pressure_bc_val = " << d_pressure_bc_val << "\n\n";

  out << "# Tissue pressure initial condition value \n";
  out << "# Default: 0\n";
  out << "tissue_pressure_ic_val = " << d_pressure_ic_val << "\n\n";

  out << "# Tissue pressure bc in north? \n";
  out << "# Default: false\n";
  out << "bc_tissue_pressure_north = " << d_pressure_bc_north << "\n\n";

  out << "# Tissue pressure bc in south? \n";
  out << "# Default: false\n";
  out << "bc_tissue_pressure_south = " << d_pressure_bc_south << "\n\n";

  out << "# Tissue pressure bc in east? \n";
  out << "# Default: false\n";
  out << "bc_tissue_pressure_east = " << d_pressure_bc_east << "\n\n";

  out << "# Tissue pressure bc in west? \n";
  out << "# Default: false\n";
  out << "bc_tissue_pressure_west = " << d_pressure_bc_west << "\n\n";
}
