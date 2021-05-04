////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_INP_H
#define UTIL_INP_H

#include "libmesh/point.h"
#include <string>
#include <vector>

using namespace libMesh;

namespace util {

struct ModelDeck {
  std::string d_model_name;
  unsigned int d_dim;
  std::string d_domain_type;
  std::vector<double> d_domain_params;
  unsigned int d_assembly_method;

  std::string d_test_name;

  std::string d_scheme_name;

  bool d_advection_active;

  bool d_decouple_nutrients;

  bool d_coupled_1d3d;

  bool d_solve_ecm;

  bool d_solve_pres_with_net_update;

  int d_seed;

  explicit ModelDeck(const std::string &filename = "")
      : d_dim(2), d_domain_type("hyper_cuboid"),
        d_domain_params(std::vector<double>(6, 0.)), d_assembly_method(2),
        d_test_name(""), d_advection_active(false), d_decouple_nutrients(false), d_seed(-1), d_coupled_1d3d(false), d_solve_ecm(true),
        d_solve_pres_with_net_update(false), d_scheme_name("solve_explicit") {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct RestartDeck {

  bool d_restart;
  std::string d_mesh_restart_file;
  std::string d_sol_restart_file;

  explicit RestartDeck(const std::string &filename = "") : d_restart(false) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct MeshDeck {
  std::string d_mesh_filename;
  bool d_read_mesh_flag;
  unsigned int d_num_elems;
  double d_mesh_size;
  bool d_use_mesh_size_for_disc;
  double d_elem_face_size;
  double d_elem_size;
  double d_face_by_h;

  explicit MeshDeck(const std::string &filename = "")
      : d_read_mesh_flag(false), d_num_elems(0), d_mesh_size(0.),
        d_use_mesh_size_for_disc(false),
        d_elem_face_size(0.), d_elem_size(0.), d_face_by_h(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct TimeDeck {
  double d_dt;
  double d_init_time;
  double d_final_time;
  unsigned int d_max_time_steps;
  unsigned int d_steps;
  unsigned int d_init_step;

  explicit TimeDeck(const std::string &filename = "")
      : d_dt(0.), d_init_time(0.), d_final_time(0.), d_max_time_steps(1000),
        d_steps(0), d_init_step(0) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct OutputDeck {
  bool d_perform_output;
  bool d_restart_save;
  unsigned int d_dt_output_interval;
  unsigned int d_dt_restart_save_interval;
  int d_debug_lvl;
  bool d_quiet;
  std::string d_outfile_tag;
  std::string d_outfilename;
  std::string d_outfilename_net;
  std::string d_output_path;
  std::string d_log_path;

  explicit OutputDeck(const std::string &filename = "")
      : d_perform_output(true), d_restart_save(false), d_dt_output_interval(1),
        d_dt_restart_save_interval(1), d_debug_lvl(2), d_quiet(false) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct SolverDeck {
  unsigned int d_linear_max_iters;
  double d_linear_tol;
  unsigned int d_nonlin_max_iters;
  double d_nonlin_tol;

  bool d_project_solution_to_physical_range;

  std::vector<std::string> d_project_fields;

  explicit SolverDeck(const std::string &filename = "")
      : d_linear_max_iters(0), d_linear_tol(0.), d_nonlin_max_iters(0),
        d_nonlin_tol(0.), d_project_solution_to_physical_range(false){

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct NutrientDeck {

  double d_lambda_P;
  double d_lambda_A;
  double d_lambda_Ph;
  double d_D_sigma;
  double d_delta_sigma;
  double d_chi_c;

  std::vector<double> d_nut_source_center;
  double d_nut_source_radius;

  explicit NutrientDeck(const std::string &filename = "")
      : d_lambda_P(0.), d_lambda_A(0.), d_lambda_Ph(0.), d_D_sigma(0.),
        d_delta_sigma(0.), d_chi_c(0.),
        d_nut_source_center(std::vector<double>(3, 0.)),
        d_nut_source_radius(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct TumorDeck {

  double d_bar_M_P;
  double d_bar_E_phi_T;
  double d_bar_E_phi_P;
  double d_epsilon_T;
  double d_epsilon_P;

  explicit TumorDeck(const std::string &filename = "")
      : d_bar_M_P(0.), d_bar_E_phi_T(0.), d_epsilon_T(0.), d_bar_E_phi_P(0.), d_epsilon_P(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct ProlificDeck {
  unsigned int d_pro_noise_num_eigenfunctions;
  unsigned int d_pro_noise_seed;
  double d_pro_noise_scale;
  double d_pro_noise_lower_bound;
  double d_pro_noise_upper_bound;
  bool d_pro_substract_avg_stoch;

  explicit ProlificDeck(const std::string &filename = "")
      : d_pro_noise_num_eigenfunctions(0),
        d_pro_noise_seed(4242),
        d_pro_noise_scale(0),
        d_pro_noise_lower_bound(0),
        d_pro_noise_upper_bound(1), d_pro_substract_avg_stoch(false) {
    if (!filename.empty())
      read_parameters(filename);
  }

  void read_parameters(const std::string &filename);

  void print(unsigned int level = 0);
};

struct HypoxicDeck {
  double d_bar_M_H;
  double d_lambda_HP;
  double d_lambda_PH;
  double d_lambda_HN;
  double d_sigma_PH;
  double d_sigma_HP;
  double d_sigma_HN;

  double d_bar_E_phi_H;
  double d_epsilon_H;

  unsigned int d_hyp_noise_num_eigenfunctions;
  unsigned int d_hyp_noise_seed;
  double d_hyp_noise_scale;
  double d_hyp_noise_lower_bound;
  double d_hyp_noise_upper_bound;
  bool d_hyp_substract_avg_stoch;

  explicit HypoxicDeck(const std::string &filename = "")
      : d_bar_M_H(0.), d_lambda_HP(0.), d_lambda_PH(0.), d_lambda_HN(0.),
        d_sigma_PH(0.), d_sigma_HP(0.), d_sigma_HN(0.),
        d_bar_E_phi_H(0.), d_epsilon_H(0.),
        d_hyp_noise_num_eigenfunctions(0),
        d_hyp_noise_seed(4242),
        d_hyp_noise_scale(0),
        d_hyp_noise_lower_bound(0), d_hyp_noise_upper_bound(1), d_hyp_substract_avg_stoch(false) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct NecroticDeck {

  double d_bar_M_N;

  explicit NecroticDeck(const std::string &filename = "") : d_bar_M_N(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct TAFDeck {

  double d_D_TAF;
  double d_delta_TAF;
  double d_lambda_TAF;

  // lower threshold of hypoxic concentration necessary to trigger TAF production
  double d_sigma_HTAF;

  // natural decay of the TAF production
  double d_lambda_TAF_deg;

  std::vector<int> d_taf_source_type;
  std::vector<std::vector<double>> d_taf_source_center;
  std::vector<double> d_taf_source_radius;

  explicit TAFDeck(const std::string &filename = "")
      : d_D_TAF(0.), d_delta_TAF(0.), d_lambda_TAF(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct ECMICData {

  /*! @brief IC type for ecm
   *
   * Choices are
   *
   * - spherical (4 params)
   * - elliptical (6 params)
   * - box (6 params)
   */
  std::string d_type;
  double d_val;
  std::vector<double> d_geom_params;

  ECMICData() : d_val(0.){};
};

struct ECMDeck {

  double d_lambda_ECM_D;
  double d_lambda_ECM_P;
  double d_bar_phi_ECM_P;
  double d_chi_h;

  ECMICData d_ecm_ic_data;

  explicit ECMDeck(const std::string &filename = "")
      : d_lambda_ECM_D(0.), d_lambda_ECM_P(0.), d_bar_phi_ECM_P(0.),
        d_chi_h(0.), d_ecm_ic_data(ECMICData()) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct MDEDeck {

  double d_D_MDE;
  double d_delta_MDE;
  double d_lambda_MDE_D;
  double d_lambda_MDE_P;

  double d_mde_ic_val;

  explicit MDEDeck(const std::string &filename = "")
      : d_D_MDE(0.), d_delta_MDE(0.), d_lambda_MDE_D(0.), d_lambda_MDE_P(0.),
        d_mde_ic_val(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct NutrientICDeck {

  /*! @brief IC type for nutrient */
  double d_nut_ic_value;

  explicit NutrientICDeck(const std::string &filename = "")
      : d_nut_ic_value(0) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct TumorICData {
  /*! @brief IC type for tumor
   *
   * Choices are
   *
   * - tumor_spherical
   * - tumor_elliptical
   * - tumor_hypoxic_spherical
   * - tumor_hypoxic_elliptical
   * - tumor_spherical_sharp
   */
  std::string d_ic_type;
  std::vector<double> d_ic_center;
  std::vector<double> d_tum_ic_radius;
  std::vector<double> d_hyp_ic_radius;

  TumorICData()
      : d_ic_center({0., 0., 0.}), d_tum_ic_radius({0., 0., 0.}),
        d_hyp_ic_radius({0., 0., 0.}){};

  TumorICData(const std::string &ic_type, const std::vector<double> &ic_center,
              const std::vector<double> &tum_ic_radius,
              const std::vector<double> &hyp_ic_radius)
      : d_ic_type(ic_type), d_ic_center(ic_center),
        d_tum_ic_radius(tum_ic_radius), d_hyp_ic_radius(hyp_ic_radius){};
};

struct TumorICDeck {

  std::vector<TumorICData> d_tum_ic_data;

  explicit TumorICDeck(const std::string &filename = "") {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct NutrientBCDeck {
  bool d_nutrient_bc_north;
  bool d_nutrient_bc_south;
  bool d_nutrient_bc_east;
  bool d_nutrient_bc_west;

  explicit NutrientBCDeck(const std::string &filename = "")
      : d_nutrient_bc_north(false), d_nutrient_bc_south(false),
        d_nutrient_bc_east(false), d_nutrient_bc_west(false) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct NetworkDeck {

  bool network_active;
  std::string d_network_init_file;
  int d_network_init_refinement;

  int d_num_points_length;
  int d_num_points_angle;
  double d_coupling_method_theta;
  bool d_compute_elem_weights;
  double d_assembly_factor_p_t;
  double d_assembly_factor_c_t;
  double d_identify_vein_pres;

  double d_identify_artery_radius;

  // 0 - surface coupling
  // 1 - line coupling
  // 2 - surface coupling for radius >= h/2, line coupling otherwise
  int d_coupling_3d1d_integration_method;
  bool d_disable_remove_redundant_vessel;
  double d_min_length_for_sprouting;

  // growth related params
  bool d_network_update;
  unsigned int d_network_update_interval;
  double d_network_update_taf_threshold;
  double d_log_normal_mean;
  double d_log_normal_std_dev;
  double d_net_radius_exponent_gamma;
  double d_network_bifurcate_prob;
  double d_min_radius;
  double d_sprouting_prob;
  double d_network_update_absolute_upper_threshold_1d;
  double d_network_update_absolute_upper_threshold_3d;
  double d_network_update_relative_upper_threshold_1d;
  double d_network_update_relative_upper_threshold_3d;

  bool d_extrapolate_nutrients_at_tips;

  double d_k_WSS;
  double d_k_s;
  double d_offset_tau;

  bool d_pressure_initial_guess_95_percent;

  bool d_remove_old_sprouters;

  // parameters below are not used currently
  int d_no_branch_dist;
  double d_new_vessel_max_angle;
  double d_branch_angle;
  int d_vessel_no_taf_effect_dist;
  unsigned int d_nonlocal_direction_search_num_points;
  double d_nonlocal_direction_search_length;
  double d_net_length_R_factor;
  double d_net_direction_lambda_g;
  bool d_network_local_search;
  double d_network_no_new_node_search_factor;

  explicit NetworkDeck(const std::string &filename = "")
      : network_active(false), d_net_direction_lambda_g(0.),
        d_net_length_R_factor(0.),
        d_network_update(true),
        d_network_update_interval(1),
        d_log_normal_mean(0.), d_log_normal_std_dev(0.),
        d_net_radius_exponent_gamma(1.), d_no_branch_dist(1),
        d_new_vessel_max_angle(0.), d_branch_angle(0.),
        d_network_update_taf_threshold(0.), d_vessel_no_taf_effect_dist(0),
        d_nonlocal_direction_search_num_points(0),
        d_nonlocal_direction_search_length(0.), d_network_local_search(false),
        d_network_no_new_node_search_factor(0.), d_num_points_length(2),
        d_num_points_angle(2), d_coupling_method_theta(0.5),
        d_assembly_factor_p_t(1.), d_assembly_factor_c_t(1.),
        d_identify_vein_pres(0.), d_compute_elem_weights(false),
        d_network_bifurcate_prob(0.9), d_min_radius(8.5e-3), d_sprouting_prob(0.9),
        d_extrapolate_nutrients_at_tips(false),
        d_pressure_initial_guess_95_percent(true),
        d_remove_old_sprouters(false),
        d_network_update_absolute_upper_threshold_1d(std::numeric_limits<double>::max()),
        d_network_update_absolute_upper_threshold_3d(std::numeric_limits<double>::max()),
        d_network_update_relative_upper_threshold_1d(std::numeric_limits<double>::max()),
        d_network_update_relative_upper_threshold_3d(std::numeric_limits<double>::max()),
        d_identify_artery_radius(0.), d_coupling_3d1d_integration_method(0), d_disable_remove_redundant_vessel(false),
        d_min_length_for_sprouting(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct Flow1DDeck {

  double d_init_vessel_mu;
  double d_in_pressure;
  double d_in_nutrient;
  double d_in_nutrient_vein;
  double d_blood_density;
  double d_D_sigma_v;

  double d_osmotic_sigma;

  std::string d_scenario;

  bool d_outlet_apply_neumann;
  bool d_inlet_apply_neumann;
  double d_outlet_neumann_val;
  double d_inlet_neumann_val;

  explicit Flow1DDeck(const std::string &filename = "")
      : d_init_vessel_mu(0.), d_in_pressure(0.), d_in_nutrient(0.),
        d_blood_density(1.),
        d_D_sigma_v(1.), d_in_nutrient_vein(0.), d_osmotic_sigma(0.),
        d_outlet_apply_neumann(false), d_outlet_neumann_val(0.),
        d_inlet_apply_neumann(false), d_inlet_neumann_val(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct FlowDeck {

  double d_tissue_flow_mu;
  double d_tissue_flow_K;
  double d_tissue_flow_coeff;
  double d_tissue_flow_rho;
  double d_tissue_flow_L_p;
  double d_tissue_nut_L_s;

  bool d_pressure_bc_north;
  bool d_pressure_bc_south;
  bool d_pressure_bc_east;
  bool d_pressure_bc_west;

  double d_pressure_bc_val;

  double d_pressure_ic_val;

  double d_mmhgFactor;

  double d_omega;

  int d_N_newton;

  explicit FlowDeck(const std::string &filename = "")
      : d_tissue_flow_mu(0.), d_tissue_flow_K(0.), d_tissue_flow_coeff(0.),
        d_tissue_flow_rho(1.),
        d_tissue_flow_L_p(0.), d_tissue_nut_L_s(0.),
        d_pressure_bc_north(false),
        d_pressure_bc_south(false), d_pressure_bc_east(false),
        d_pressure_bc_west(false), d_pressure_bc_val(0.),
        d_pressure_ic_val(0.), d_mmhgFactor(133.322), d_omega(0.0), d_N_newton(0) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
  void print(unsigned int level = 0);
};

struct InputDeck : public ModelDeck,
                   public RestartDeck,
                   public MeshDeck,
                   public TimeDeck,
                   public SolverDeck,
                   public OutputDeck,
                   public NutrientDeck,
                   public TumorDeck,
                   public HypoxicDeck,
                   public ProlificDeck,
                   public NecroticDeck,
                   public TAFDeck,
                   public ECMDeck,
                   public MDEDeck,
                   public NutrientICDeck,
                   public TumorICDeck,
                   public NutrientBCDeck,
                   public NetworkDeck,
                   public Flow1DDeck,
                   public FlowDeck {

public:
  explicit InputDeck(const std::string &filename = "")
      : ModelDeck(filename), RestartDeck(filename), MeshDeck(filename),
        TimeDeck(filename), SolverDeck(filename), OutputDeck(filename),
        NutrientDeck(filename), TumorDeck(filename), HypoxicDeck(filename),
        ProlificDeck(filename), NecroticDeck(filename), TAFDeck(filename),
        ECMDeck(filename), MDEDeck(filename), NutrientICDeck(filename),
        TumorICDeck(filename), NutrientBCDeck(filename), NetworkDeck(filename),
        Flow1DDeck(filename), FlowDeck(filename){};

  //
  void print(unsigned int level = 0){

    {ModelDeck *deck = this;
  deck->print();
}

{
  RestartDeck *deck = this;
  deck->print();
}

{
  MeshDeck *deck = this;
  deck->print();
}

{
  TimeDeck *deck = this;
  deck->print();
}

{
  SolverDeck *deck = this;
  deck->print();
}

{
  OutputDeck *deck = this;
  deck->print();
}

{
  NutrientDeck *deck = this;
  deck->print();
}

{
  TumorDeck *deck = this;
  deck->print();
}

{
  HypoxicDeck *deck = this;
  deck->print();
}

{
  ProlificDeck *deck = this;
  deck->print();
}

{
  NecroticDeck *deck = this;
  deck->print();
}

{
  TAFDeck *deck = this;
  deck->print();
}

{
  ECMDeck *deck = this;
  deck->print();
}

{
  MDEDeck *deck = this;
  deck->print();
}

{
  NutrientICDeck *deck = this;
  deck->print();
}

{
  TumorICDeck *deck = this;
  deck->print();
}

{
  NutrientBCDeck *deck = this;
  deck->print();
}

{
  NetworkDeck *deck = this;
  deck->print();
}

{
  Flow1DDeck *deck = this;
  deck->print();
}

{
  FlowDeck *deck = this;
  deck->print();
}
}; // namespace util
}
;

} // namespace util

#endif // UTIL_INP_H
