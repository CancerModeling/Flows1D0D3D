////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NET_INP_H
#define NET_INP_H

#include "libmesh/point.h"
#include <string>
#include <vector>

using namespace libMesh;

namespace util {

struct ModelDeck {
  unsigned int d_dim;
  std::string d_domain_type;
  std::vector<double> d_domain_params;
  unsigned int d_assembly_method;

  std::string d_test_name;

  bool d_advection_active;

  bool d_decouple_nutrients;

  explicit ModelDeck(const std::string &filename = "")
      : d_dim(2), d_domain_type("hyper_cuboid"),
        d_domain_params(std::vector<double>(6, 0.)), d_assembly_method(2),
        d_test_name(""), d_advection_active(false), d_decouple_nutrients
        (false) {

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

  explicit SolverDeck(const std::string &filename = "")
      : d_linear_max_iters(0), d_linear_tol(0.), d_nonlin_max_iters(0),
        d_nonlin_tol(0.), d_project_solution_to_physical_range(false) {

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
  double d_epsilon_T;

  explicit TumorDeck(const std::string &filename = "")
      : d_bar_M_P(0.), d_bar_E_phi_T(0.), d_epsilon_T(0.) {

    if (!filename.empty())
      read_parameters(filename);
  };

  //
  void read_parameters(const std::string &filename);

  //
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

  explicit HypoxicDeck(const std::string &filename = "")
      : d_bar_M_H(0.), d_lambda_HP(0.), d_lambda_PH(0.), d_lambda_HN(0.),
        d_sigma_PH(0.), d_sigma_HP(0.), d_sigma_HN(0.) {

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

  ECMICData() : d_val(0.) {};
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



  unsigned int d_network_update_interval;

  double d_network_update_taf_threshold;

  int d_no_branch_dist;
  double d_new_vessel_max_angle;
  double d_branch_angle;
  int d_vessel_no_taf_effect_dist;
  unsigned int d_nonlocal_direction_search_num_points;
  double d_nonlocal_direction_search_length;

  double d_log_normal_mean;
  double d_log_normal_std_dev;
  double d_net_length_R_factor;
  double d_net_radius_exponent_gamma;
  double d_net_direction_lambda_g;

  bool d_network_local_search;

  double d_network_no_new_node_search_factor;

  double d_network_bifurcate_prob;

  explicit NetworkDeck(const std::string &filename = "")
      : network_active(false), d_net_direction_lambda_g(0.),
        d_net_length_R_factor(0.), d_network_update_interval(1),
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
        d_network_bifurcate_prob(0.6) {

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

  explicit Flow1DDeck(const std::string &filename = "")
      : d_init_vessel_mu(0.), d_in_pressure(0.), d_in_nutrient(0.),
          d_blood_density(1.),
      d_D_sigma_v(1.), d_in_nutrient_vein(0.), d_osmotic_sigma(0.) {

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

  explicit FlowDeck(const std::string &filename = "")
      : d_tissue_flow_mu(0.), d_tissue_flow_K(0.), d_tissue_flow_coeff(0.),
      d_tissue_flow_rho(1.),
      d_tissue_flow_L_p(0.), d_tissue_nut_L_s(0.),
      d_pressure_bc_north(false),
        d_pressure_bc_south(false), d_pressure_bc_east(false),
        d_pressure_bc_west(false), d_pressure_bc_val(0.),
        d_pressure_ic_val(0.), d_mmhgFactor(133.322) {

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
        NecroticDeck(filename), TAFDeck(filename), ECMDeck(filename),
        MDEDeck(filename), NutrientICDeck(filename), TumorICDeck(filename),
        NutrientBCDeck(filename), NetworkDeck(filename), Flow1DDeck(filename),
        FlowDeck(filename){};

  //
  void print(unsigned int level = 0) {

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
};
};

} // namespace util

#endif // NET_INP_H