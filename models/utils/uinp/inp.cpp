////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "inp.hpp"
#include "csv.hpp"
#include "libmesh/getpot.h"

void util::ModelDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_model_name = input("model_name", "");
  if (d_model_name.empty())
    libmesh_error_msg("Model name is not valid");

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

  d_assembly_method = input("assembly_method", 2);

  d_test_name = input("test_name", "");

  d_scheme_name = input("scheme_name", "");

  d_advection_active = input("advection_active", true);

  d_decouple_nutrients = input("network_decouple_nutrients", true);

  d_seed = input("seed", -1);

  d_coupled_1d3d = input("coupled_1d3d", false);

  d_solve_ecm = input("solve_ecm", true);

  d_solve_pres_with_net_update = input("solve_pres_with_net_update", false);
}

void util::ModelDeck::print(unsigned int level) {}

void util::RestartDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_restart = input("restart", false);
  d_mesh_restart_file = input("mesh_restart_file", "");
  d_sol_restart_file = input("solution_restart_file", "");
}

void util::RestartDeck::print(unsigned int level) {}

void util::MeshDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_mesh_filename = input("mesh_name", "");
  d_num_elems = input("mesh_n_elements", 0);
  if (d_num_elems == 0) {
    // check if mesh size is provided
    d_mesh_size = input("mesh_size", 0.);

    if (d_mesh_size < 1.0e-12) {
      libmesh_error_msg("Number of elements for 2D/3D mesh can not be zero.");
      exit(1);
    } else {
      d_use_mesh_size_for_disc = true;
    }
  }
  d_read_mesh_flag = input("mesh_read_file", false);
}

void util::MeshDeck::print(unsigned int level) {}

void util::TimeDeck::read_parameters(const std::string &filename) {

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

void util::TimeDeck::print(unsigned int level) {}

void util::OutputDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_perform_output = input("perform_output", true);
  d_dt_output_interval = input("output_interval", 1);

  d_restart_save = input("restart_save", false);
  d_dt_restart_save_interval = input("restart_save_interval", 1);

  bool dbg_out = input("output_debug_info", true);
  d_quiet = !dbg_out;

  d_output_path = input("output_path", "");
  d_log_path = input("log_output_path", "");
  d_outfile_tag = input("output_tag", "");

  d_outfilename = d_output_path + "output";
  d_outfilename_net = d_output_path + "net_output";

  // For now do not put tag on output files so that paraview state
  // could be reused for any simulation
  //  if (!d_outfile_tag.empty()) {
  //    d_outfilename += "_" + d_outfile_tag;
  //    d_outfilename_net += "_" + d_outfile_tag;
  //  }
}

void util::OutputDeck::print(unsigned int level) {}

void util::SolverDeck::read_parameters(const std::string &filename) {

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

  // find which fields need to be projected
  bool field_project = input("project_prolific", false);
  if (field_project)
    d_project_fields.emplace_back("prolific");

  field_project = input("project_hypoxic", false);
  if (field_project)
    d_project_fields.emplace_back("hypoxic");

  field_project = input("project_necrotic", false);
  if (field_project)
    d_project_fields.emplace_back("necrotic");

  field_project = input("project_tumor", false);
  if (field_project)
    d_project_fields.emplace_back("tumor");
}

void util::SolverDeck::print(unsigned int level) {}

void util::NutrientDeck::read_parameters(const std::string &filename) {

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

  d_nut_source_center[0] = input("nut_source_center_x", 0.);
  d_nut_source_center[1] = input("nut_source_center_y", 0.);
  d_nut_source_center[2] = input("nut_source_center_z", 0.);
  d_nut_source_radius = input("nut_source_radius", 0.);
}

void util::NutrientDeck::print(unsigned int level) {}

void util::TumorDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_bar_M_P = input("bar_M_P", 200.);
  d_bar_E_phi_T = input("bar_E_phi_T", 0.045);
  d_epsilon_T = input("epsilon_T", 0.005);

  d_bar_E_phi_P = input("bar_E_phi_P", 0.0);
  d_epsilon_P = input("epsilon_P", 0.005);
}

void util::TumorDeck::print(unsigned int level) {}

void util::HypoxicDeck::read_parameters(const std::string &filename) {

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

  d_bar_E_phi_H = input("bar_E_phi_H", 0.0);
  d_epsilon_H = input("epsilon_H", 0.005);

  d_hyp_noise_num_eigenfunctions = input("hyp_noise_num_eigenfunctions", 0);
  d_hyp_noise_seed = input("hyp_noise_seed", 104020);
  d_hyp_noise_scale = input("hyp_noise_scale", 0.1);
  d_hyp_noise_lower_bound = input("hyp_noise_lower_bound", 0.0);
  d_hyp_noise_upper_bound = input("hyp_noise_upper_bound", 1.0);
  d_hyp_substract_avg_stoch = input("hyp_substract_avg_stoch", false);
}

void util::HypoxicDeck::print(unsigned int level) {}

void util::ProlificDeck::read_parameters(const std::string &filename) {
  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);
  d_pro_noise_num_eigenfunctions = input("pro_noise_num_eigenfunctions", 0);
  d_pro_noise_seed = input("pro_noise_seed", 104020);
  d_pro_noise_scale = input("pro_noise_scale", 0.1);
  d_pro_noise_lower_bound = input("pro_noise_lower_bound", 0.0);
  d_pro_noise_upper_bound = input("pro_noise_upper_bound", 1.0);
  d_pro_substract_avg_stoch = input("pro_substract_avg_stoch", false);
}

void util::ProlificDeck::print(unsigned int level) {}

void util::NecroticDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_bar_M_N = input("bar_M_N", 0.);
}

void util::NecroticDeck::print(unsigned int level) {}

void util::TAFDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_D_TAF = input("D_TAF", 10.);
  d_delta_TAF = input("delta_TAF", 1.);
  d_lambda_TAF = input("lambda_TAF", 10.);
  d_sigma_HTAF = input("sigma_HTAF", 0.);
  d_lambda_TAF_deg = input("lambda_TAF_deg", 0.);

  // see if taf source file is provided
  bool read_csv = false;
  std::string csv_file = input("taf_source_file", "");
  if (!csv_file.empty())
    read_csv = true;

  if (read_csv) {

    // header
    // type (int), center (x,y,z), tum_radius(r1, r2, r3), hyp_radius(r1, r2,
    // r3)
    io::CSVReader<5> in(csv_file);

    in.read_header(io::ignore_extra_column, "type", "cx", "cy", "cz", "r");

    int type;
    double cx, cy, cz, r;
    while (in.read_row(type, cx, cy, cz, r)) {
      d_taf_source_type.push_back(type);
      d_taf_source_center.push_back({cx, cy, cz});
      d_taf_source_radius.push_back(r);
    }
  }
}

void util::TAFDeck::print(unsigned int level) {}

void util::ECMDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_lambda_ECM_D = input("lambda_ECM_D", 1.);
  d_lambda_ECM_P = input("lambda_ECM_P", 1.);
  d_bar_phi_ECM_P = input("bar_phi_ECM_P", 0.5);
  d_chi_h = input("chi_h", 0.);

  // read ic information for ecm

  // locate the domain of stroma
  d_ecm_ic_data.d_type = input("ECM_ic_domain_type", "");

  // get the initial value of ecm in stroma
  d_ecm_ic_data.d_val = input("ECM_ic_val", 1.);

  // get the geometrical parameters
  unsigned int num_params = input("ECM_ic_num_params", 0);

  for (unsigned int i = 0; i < num_params; i++)
    d_ecm_ic_data.d_geom_params.push_back(
      input("ECM_ic_params_" + std::to_string(i + 1), 0.));
}

void util::ECMDeck::print(unsigned int level) {}

void util::MDEDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_D_MDE = input("D_MDE", 10.);
  d_delta_MDE = input("delta_MDE", 1.);
  d_lambda_MDE_D = input("lambda_MDE_D", 1.);
  d_lambda_MDE_P = input("lambda_MDE_P", 1.);

  d_mde_ic_val = input("MDE_ic_val", 0.);
}

void util::MDEDeck::print(unsigned int level) {}

void util::NutrientICDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_nut_ic_value = input("ic_nutrient_value", 0.5);
  // out << "nutrient ic = " << d_nut_ic_value;
}

void util::NutrientICDeck::print(unsigned int level) {}

void util::TumorICDeck::read_parameters(const std::string &filename) {

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

    in.read_header(io::ignore_extra_column, "type", "cx", "cy", "cz", "tum_rx",
                   "tum_ry", "tum_rz", "hyp_rx", "hyp_ry", "hyp_rz");

    int type;
    double cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz;
    while (in.read_row(type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry,
                       hyp_rz)) {

      std::string ic_type;
      if (type == 1)
        ic_type = "tumor_spherical";
      else if (type == 2)
        ic_type = "tumor_elliptical";
      else if (type == 3)
        ic_type = "tumor_hypoxic_spherical";
      else if (type == 4)
        ic_type = "tumor_hypoxic_elliptical";
      else if (type == 5)
        ic_type = "tumor_spherical_sharp";

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

      auto type = input("ic_tumor_type" + postfix, 0);
      if (type == 1)
        data.d_ic_type = "tumor_spherical";
      else if (type == 2)
        data.d_ic_type = "tumor_elliptical";
      else if (type == 3)
        data.d_ic_type = "tumor_hypoxic_spherical";
      else if (type == 4)
        data.d_ic_type = "tumor_hypoxic_elliptical";
      else if (type == 5)
        data.d_ic_type = "tumor_spherical_sharp";

      data.d_ic_center[0] = input("ic_tumor_center_x" + postfix, 0.);
      data.d_ic_center[1] = input("ic_tumor_center_y" + postfix, 0.);
      data.d_ic_center[2] = input("ic_tumor_center_z" + postfix, 0.);
      data.d_tum_ic_radius[0] = input("ic_tumor_radius_x" + postfix, 0.);
      data.d_tum_ic_radius[1] = input("ic_tumor_radius_y" + postfix, 0.);
      data.d_tum_ic_radius[2] = input("ic_tumor_radius_z" + postfix, 0.);
      data.d_hyp_ic_radius[0] = input("ic_hypoxic_radius_x" + postfix, 0.);
      data.d_hyp_ic_radius[1] = input("ic_hypoxic_radius_y" + postfix, 0.);
      data.d_hyp_ic_radius[2] = input("ic_hypoxic_radius_z" + postfix, 0.);

      // add data
      d_tum_ic_data[i] = data;
    }
  }
}

void util::TumorICDeck::print(unsigned int level) {}

void util::NutrientBCDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  d_nutrient_bc_north = input("bc_nutrient_north", false);
  d_nutrient_bc_south = input("bc_nutrient_south", false);
  d_nutrient_bc_east = input("bc_nutrient_east", false);
  d_nutrient_bc_west = input("bc_nutrient_west", false);
}

void util::NutrientBCDeck::print(unsigned int level) {}

void util::NetworkDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  // get model name first
  std::string model_name = input("model_name", "");
  if (model_name == "TwoSpecies" or model_name == "AvaFV")
    return;

  network_active = input("is_network_active", false);
  d_network_init_file = input("network_init_file", "");
  d_network_init_refinement = input("network_init_refinement", 1);
  d_num_points_length = input("network_discret_cyl_length", 2);
  d_num_points_angle = input("network_discret_cyl_angle", 2);
  d_coupling_method_theta = input("network_coupling_method_theta", 0.5);
  d_compute_elem_weights = input("network_compute_elem_weights", true);

  if (input.have_variable("assembly_factor_p_t"))
    d_assembly_factor_p_t = input("assembly_factor_p_t", 1.);
  else if (input.have_variable("coupling_factor_p_t"))
    d_assembly_factor_p_t = input("coupling_factor_p_t", 1.);

  d_assembly_factor_c_t = input("assembly_factor_c_t", 1.);
  d_identify_vein_pres = input("identify_vein_pressure", 0.);
  d_identify_artery_radius = input("identify_artery_radius", -1.);
  if (d_identify_artery_radius < 0.)
    libmesh_error_msg("Error. Must specify radius to identify artery.");

  d_extrapolate_nutrients_at_tips = input("extrapolate_nutrients_at_tips", true);

  d_coupling_3d1d_integration_method = input("coupling_3d1d_integration_method", 0);
  d_disable_remove_redundant_vessel = input("disable_remove_redundant_vessel", false);
  d_min_length_for_sprouting = input("min_length_for_sprouting", 0.);

  // growth related params
  d_network_update = input("network_update", true);
  d_network_update_interval = input("network_update_interval", 3);
  //   if (d_network_update and d_network_update_interval != 3) {
  //     d_network_update_interval = 3;
  //     libmesh_warning("Currently network update interval is fixed to 3 and can not be changed from input file.");
  //   }

  d_network_update_taf_threshold = input("network_update_taf_threshold", 0.);
  d_log_normal_mean = input("log_normal_mean", 0.);
  d_log_normal_std_dev = input("log_normal_std_dev", 0.);
  d_net_radius_exponent_gamma = input("network_radius_exponent_gamma", 1.);
  d_network_bifurcate_prob = input("network_bifurcate_probability", 0.9);
  d_min_radius = input("network_min_radius", 8.5e-3);
  d_sprouting_prob = input("network_sprouting_prob", 0.9);

  d_network_update_absolute_upper_threshold_1d = input("network_update_absolute_upper_threshold_1d", std::numeric_limits<double>::max());
  d_network_update_absolute_upper_threshold_3d = input("network_update_absolute_upper_threshold_3d", std::numeric_limits<double>::max());
  d_network_update_relative_upper_threshold_1d = input("network_update_relative_upper_threshold_1d", std::numeric_limits<double>::max());
  d_network_update_relative_upper_threshold_3d = input("network_update_relative_upper_threshold_3d", std::numeric_limits<double>::max());

  // parameters which are not used currently
  d_no_branch_dist = input("network_no_branch_dist", 1);
  d_new_vessel_max_angle = input("network_new_veesel_max_angle", M_PI / 4.);
  d_branch_angle = input("network_branch_angle", M_PI / 8.);
  d_vessel_no_taf_effect_dist = input("network_vessel_no_taf_dist", 5);
  d_nonlocal_direction_search_num_points =
    input("network_nonlocal_search_num_points", 2);
  d_nonlocal_direction_search_length =
    input("network_nonlocal_search_length_factor", 10.);
  d_network_local_search = input("network_local_search", true);
  d_network_no_new_node_search_factor =
    input("network_no_new_node_search_factor", 0.5);
  d_net_direction_lambda_g = input("vessel_lambda_g", 0.);
  d_net_length_R_factor = input("vessel_R_factor", 0.);

  d_k_WSS = input("k_WSS", 0.45);
  d_k_s = input("k_s", 0.25);
  d_offset_tau = input("offset_tau", 0.02);

  d_remove_old_sprouters = input("remove_old_sprouters", false);

  d_pressure_initial_guess_95_percent = input("pressure_initial_guess_95_percent", true);
}

void util::NetworkDeck::print(unsigned int level) {}

void util::Flow1DDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  // get model name first
  std::string model_name = input("model_name", "");
  if (model_name == "TwoSpecies" or model_name == "AvaFV")
    return;

  d_init_vessel_mu = input("init_vessel_viscosity", 1.);
  d_in_pressure = input("vessel_in_pressure", 0.);
  d_in_nutrient = input("vessel_in_nutrient", 0.);
  d_in_nutrient_vein = input("vessel_in_nutrient_vein", 0.);
  d_blood_density = input("vessel_blood_density", 1.);
  d_D_sigma_v = input("vessel_D_sigma", 1.);

  d_osmotic_sigma = input("osmotic_reflection_coeff", 0.1);

  d_scenario = input("scenario", "scenario not available");

  d_outlet_apply_neumann = input("outlet_apply_neumann", false);
  d_inlet_apply_neumann = input("inlet_apply_neumann", false);
  d_outlet_neumann_val = input("outlet_neumann_val", 0.);
  d_inlet_neumann_val = input("inlet_neumann_val", 0.);

  // we can not have neumann on both inlet and outlet
  if (d_outlet_apply_neumann and d_inlet_apply_neumann)
    libmesh_error_msg("Error: Can not apply Neumann boundary condition on "
                      "both inlets and outlets.");
}

void util::Flow1DDeck::print(unsigned int level) {}

void util::FlowDeck::read_parameters(const std::string &filename) {

  // Open file with model setup
  if (filename.empty())
    return;

  GetPot input(filename);

  // get model name first
  std::string model_name = input("model_name", "");
  if (model_name == "TwoSpecies" or model_name == "AvaFV")
    return;

  d_tissue_flow_mu = input("tissue_flow_viscosity", 1.);
  d_tissue_flow_K = input("tissue_flow_K", 1.);
  d_tissue_flow_coeff = d_tissue_flow_K / d_tissue_flow_mu;
  d_tissue_flow_rho = input("tissue_flow_density", 1.);
  d_tissue_flow_L_p = input("tissue_flow_L_p", 1.);
  d_tissue_nut_L_s = input("tissue_nut_L_s", 1.);

  d_pressure_bc_val = input("tissue_pressure_bc_val", 0.);
  d_pressure_ic_val = input("tissue_pressure_ic_val", 0.);

  d_pressure_bc_north = input("bc_tissue_pressure_north", false);
  d_pressure_bc_south = input("bc_tissue_pressure_south", false);
  d_pressure_bc_east = input("bc_tissue_pressure_east", false);
  d_pressure_bc_west = input("bc_tissue_pressure_west", false);

  d_mmhgFactor = input("mmhg_factor", 1.);

  d_N_newton = input("N_newton", 0);
  d_omega = input("omega", 0.0);
}

void util::FlowDeck::print(unsigned int level) {}
