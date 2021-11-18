////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_3d_solver.hpp"
#include "heart_to_breast_1d_solver.hpp"
#include "random_dist.hpp"
#include "tree_search.hpp"
#include "vtk_io_libmesh.hpp"
#include "vtk_writer.hpp"
#include "libmesh_utils.hpp"
#include "tree_search.hpp"
#include <cfloat>
#include <fmt/format.h>

namespace macrocirculation {

namespace {
// creates outlet at random location and assigns outlet radii randomly
void set_perfusion_pts(std::string out_dir,
                       int num_pts,
                       std::vector<lm::Point> &pts,
                       std::vector<double> &radii,
                       lm::EquationSystems &eq_sys,
                       HeartToBreast3DSolver &model) {

  // initialize random number generator
  int seed = 0;
  srand(seed);

  const auto &input = eq_sys.parameters.get<HeartToBreast3DSolverInputDeck *>("input_deck");
  const auto &mesh = eq_sys.get_mesh();

  // get length of the domain
  auto bbox = lm::MeshTools::create_bounding_box(mesh);
  auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
  auto l = (bbox.min() - bbox.max()).norm();

  // create list of element centers for tree search
  int nelems = mesh.n_elem();
  std::vector<lm::Point> elem_centers(nelems, lm::Point());
  for (const auto &elem : mesh.element_ptr_range())
    elem_centers[elem->id()] = elem->centroid();

  // randomly select desired number of element centers as outlet perfusion points
  int npts = num_pts;
  pts.resize(npts);
  radii.resize(npts);
  std::vector<int> sel_elems;
  double min_dist = l / npts;
  for (int i = 0; i < 10 * npts; i++) {
    if (sel_elems.size() == npts)
      break;

    // get random integer between 0 and nelems - 1
    int e = rand() % (nelems - 1);
    if (locate_in_set(e, sel_elems) != -1)
      continue;

    // e is not in existing list so check if it is a good candidate
    bool not_suitable = false;
    for (auto ee : sel_elems) {
      auto dx = elem_centers[ee] - elem_centers[e];
      if (dx.norm() < min_dist) {
        not_suitable = true;
        break;
      }
    }

    if (not_suitable)
      continue;

    // add element to the list
    sel_elems.push_back(e);
  }

  // if at this point we do not have enough elements in sel_elems than exit
  if (sel_elems.size() < npts) {
    std::cerr << "Error: Increase threshold for creating random points for perfusion\n";
    exit(EXIT_FAILURE);
  }

  // add cooardinates and radius (based on uniform distribution)
  DistributionSample<UniformDistribution> uni_dist(min_dist / 10., min_dist / 3., seed);
  for (int i = 0; i < npts; i++) {
    pts[i] = elem_centers[sel_elems[i]];
    radii[i] = uni_dist();
  }
}
} // namespace

// input class definitions

HeartToBreast3DSolverInputDeck::HeartToBreast3DSolverInputDeck(const std::string &filename)
    : d_rho_cap(1.), d_rho_tis(1.), d_K_cap(1.e-5), d_K_tis(1.e-11),
      d_Lp_art_cap(1.e-6), d_Lp_vein_cap(1.e-8), d_Lp_cap_tis(1e-11),
      d_Dnut_cap(1e-3), d_Dtis_cap(1.e-6), d_Lnut_cap_tis(0.01),
      d_N_bar_cap(1e2), d_N_bar_surf_cap(1.e-2),
      d_rnut_cap(0.), d_rnut_art_cap(0.), d_rnut_vein_cap(1.),
      d_lambda_P(5.), d_lambda_A(0.), d_tum_mob(1.),
      d_tum_dw(0.45), d_tum_eps(0.0158),
      d_T(1.), d_dt(0.01), d_h(0.1), d_mesh_file(""), d_out_dir(""),
      d_perf_regularized(false),
      d_perf_fn_type("const"), d_perf_neigh_size({1., 4.}),
      d_debug_lvl(0) {
  if (!filename.empty())
    read_parameters(filename);
}

void HeartToBreast3DSolverInputDeck::read_parameters(const std::string &filename) {
  GetPot input(filename);
  d_rho_cap = input("rho_cap", d_rho_cap);
  d_rho_tis = input("rho_tis", d_rho_tis);
  d_K_cap = input("K_cap", d_K_cap);
  d_K_tis = input("K_tis", d_K_tis);
  d_Lp_art_cap = input("Lp_art_cap", d_Lp_art_cap);
  d_Lp_vein_cap = input("Lp_vein_cap", d_Lp_vein_cap);
  d_Lp_cap_tis = input("Lp_cap_tis", d_Lp_cap_tis);
  d_Dnut_cap = input("Dnut_cap", d_Dnut_cap);
  d_Dtis_cap = input("Dtis_cap", d_Dtis_cap);
  d_Lnut_cap_tis = input("Lnut_cap_tis", d_Lnut_cap_tis);
  d_N_bar_cap = input("N_bar_cap", d_N_bar_cap);
  d_N_bar_surf_cap = input("N_bar_surf_cap", d_N_bar_surf_cap);
  d_rnut_cap = input("rnut_cap", d_rnut_cap);
  d_rnut_art_cap = input("rnut_art_cap", d_rnut_art_cap);
  d_rnut_vein_cap = input("rnut_vein_cap", d_rnut_vein_cap);
  d_lambda_P = input("lambda_P", d_lambda_P);
  d_lambda_A = input("lambda_A", d_lambda_A);
  d_tum_mob = input("tum_mob", d_tum_mob);
  d_tum_dw = input("tum_dw", d_tum_dw);
  d_tum_eps = input("tum_eps", d_tum_eps);
  d_T = input("T", d_T);
  d_dt = input("dt", d_dt);
  d_h = input("h", d_h);
  d_mesh_file = input("mesh_file", d_mesh_file);
  d_out_dir = input("out_dir", d_out_dir);
  d_perf_regularized = input("regularized_source", d_perf_regularized ? 0 : 1) == 0;
  d_perf_fn_type = input("perf_fn_type", d_perf_fn_type);
  d_perf_neigh_size.first = input("perf_neigh_size_min", d_perf_neigh_size.first);
  d_perf_neigh_size.second = input("perf_neigh_size_max", d_perf_neigh_size.second);
  d_debug_lvl = input("debug_lvl", d_debug_lvl);
}

std::string HeartToBreast3DSolverInputDeck::print_str() {
  std::ostringstream oss;
  oss << "rho_cap = " << d_rho_cap << "\n";
  oss << "rho_tis = " << d_rho_tis << "\n";
  oss << "K_cap = " << d_K_cap << "\n";
  oss << "K_tis = " << d_K_tis << "\n";
  oss << "L_art_cap = " << d_Lp_art_cap << "\n";
  oss << "Lp_cap_tis = " << d_Lp_cap_tis << "\n";
  oss << "N_bar_cap = " << d_N_bar_cap << "\n";
  oss << "T = " << d_T << "\n";
  oss << "dt = " << d_dt << "\n";
  oss << "h = " << d_h << "\n";
  oss << "mesh_file = " << d_mesh_file << "\n";
  oss << "out_dir = " << d_out_dir << "\n";
  return oss.str();
}

// solver class definitions
HeartToBreast3DSolver::HeartToBreast3DSolver(MPI_Comm mpi_comm,
                                             lm::Parallel::Communicator *libmesh_comm,
                                             HeartToBreast3DSolverInputDeck &input,
                                             lm::ReplicatedMesh &mesh,
                                             lm::EquationSystems &eq_sys,
                                             lm::TransientLinearImplicitSystem &p_cap,
                                             lm::TransientLinearImplicitSystem &p_tis,
                                             lm::TransientLinearImplicitSystem &nut_cap,
                                             lm::TransientLinearImplicitSystem &nut_tis,
                                             lm::TransientLinearImplicitSystem &tum,
                                             lm::ExplicitSystem &K_tis_field,
                                             lm::ExplicitSystem &Dnut_tis_field,
                                             lm::ExplicitSystem &N_bar_cap_field,
                                             lm::ExplicitSystem &N_bar_sruf_cap_field,
                                             Logger &log)
    : BaseModel(libmesh_comm, mesh, eq_sys, log, "HeartToBreast3DSolver"),
      d_input(input),
      d_p_cap(this, d_mesh, p_cap),
      d_p_tis(this, d_mesh, p_tis),
      d_nut_cap(this, d_mesh, nut_cap),
      d_nut_tis(this, d_mesh, nut_tis),
      d_tum(this, d_mesh, tum),
      d_K_tis_field(K_tis_field),
      d_Dnut_tis_field(Dnut_tis_field),
      d_N_bar_cap_field(N_bar_cap_field),
      d_N_bar_surf_cap_field(N_bar_sruf_cap_field) {

  d_nut_convert_factor = 5.550930*1.e-5;
  d_dt = input.d_dt;
  d_log("created HeartToBreast3DSolver object\n");
}

void HeartToBreast3DSolver::write_perfusion_output(std::string out_file) {

  auto vtu_writer = VTKWriter(out_file);
  add_points(d_perf_pts, vtu_writer.d_d_p);
  add_array("Radius", d_perf_radii, vtu_writer.d_d_p);
  add_array("Ball_Radius", d_perf_ball_radii, vtu_writer.d_d_p);
  add_array("p_art_outlet", d_perf_pres, vtu_writer.d_d_p);
  add_array("p_vein_outlet", d_perf_pres_vein, vtu_writer.d_d_p);
  add_array("nut_art_outlet", d_perf_nut, vtu_writer.d_d_p);
  add_array("nut_vein_outlet", d_perf_nut_vein, vtu_writer.d_d_p);
  add_array("p_cap_outlet", d_perf_p_3d_weighted, vtu_writer.d_d_p);
  add_array("nut_cap_outlet", d_perf_nut_3d_weighted, vtu_writer.d_d_p);
  vtu_writer.write();
}

void HeartToBreast3DSolver::setup() {
  // TODO setup other things
  //d_eq_sys.init();
}
double HeartToBreast3DSolver::get_time() const {
  return d_time;
}
void HeartToBreast3DSolver::solve(bool solve_tum) {
  auto solve_clock = std::chrono::steady_clock::now();
  d_p_cap.solve();
  d_log("capillary pressure solve time = " + std::to_string(time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

  solve_clock = std::chrono::steady_clock::now();
  d_p_tis.solve();
  d_log("tissue pressure solve time = " + std::to_string(time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

  solve_clock = std::chrono::steady_clock::now();
  d_nut_cap.solve();
  d_log("capillary nutrient solve time = " + std::to_string(time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

  solve_clock = std::chrono::steady_clock::now();
  d_nut_tis.solve();
  d_log("tissue nutrient solve time = " + std::to_string(time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

  if (solve_tum) {
    solve_clock = std::chrono::steady_clock::now();
    d_tum.solve();
    d_log("tumor solve time = " + std::to_string(time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");
  }
}
void HeartToBreast3DSolver::write_output() {
  static int out_n = 0;
  // out 3d data
  VTKIO(d_mesh).write_equation_systems(d_input.d_out_dir + "/output_3D_" + std::to_string(out_n) + ".pvtu", d_eq_sys);
  // out 1d data
  write_perfusion_output(d_input.d_out_dir + "/output_3D_perf_" + std::to_string(out_n) + ".vtu");
  // out 3d qoi data
  std::vector<double> qoi;
  qoi.push_back(d_p_cap.compute_qoi("linf"));
  qoi.push_back(d_p_cap.compute_qoi("l1"));
  qoi.push_back(d_p_cap.compute_qoi("l2"));
  qoi.push_back(d_p_tis.compute_qoi("linf"));
  qoi.push_back(d_p_tis.compute_qoi("l1"));
  qoi.push_back(d_p_tis.compute_qoi("l2"));
  qoi.push_back(d_nut_cap.compute_qoi("linf"));
  qoi.push_back(d_nut_cap.compute_qoi("l1"));
  qoi.push_back(d_nut_cap.compute_qoi("l2"));
  qoi.push_back(d_nut_tis.compute_qoi("linf"));
  qoi.push_back(d_nut_tis.compute_qoi("l1"));
  qoi.push_back(d_nut_tis.compute_qoi("l2"));
  qoi.push_back(d_tum.compute_qoi("linf", 0));
  qoi.push_back(d_tum.compute_qoi("l1", 0));
  qoi.push_back(d_tum.compute_qoi("l2", 0));

  if (d_procRank == 0) {
    std::string fn = fmt::format("{}qoi_3d.txt", d_input.d_out_dir);
    if (out_n == 0) {
      std::ofstream of;
      of.open(fn);
      of << "t, p_cap_linf, p_cap_l1, p_cap_l2, p_tis_linf, p_tis_l1, p_tis_l2, nut_cap_linf, nut_cap_l1, "
            "nut_cap_l2, nut_tis_linf, nut_tis_l1, nut_tis_l2, tum_linf, tum_l1, tum_l2\n";
    }
    std::ofstream of;
    of.open(fn, std::ios_base::app);
    of << d_time << ", ";
    for (size_t i = 0; i < qoi.size(); i++)
      of << qoi[i] << (i < qoi.size() - 1 ? ", " : "\n");
    of.close();
  }

  out_n++;
}
void HeartToBreast3DSolver::set_output_folder(std::string output_dir) {
}
void HeartToBreast3DSolver::setup_1d3d(const std::vector<VesselTipCurrentCouplingData> &data_1d) {
  if (d_input.d_perf_regularized) {
    std::cout << "Setting up regularized perfusion sources\n";
    setup_1d3d_reg_source(data_1d);
  } else {
    std::cout << "Setting up uniform partitioned perfusion sources\n";
    setup_1d3d_partition(data_1d);
  }
}

void HeartToBreast3DSolver::setup_1d3d_reg_source(const std::vector<VesselTipCurrentCouplingData> &data_1d) {

  auto num_perf_outlets = data_1d.size();
  if (num_perf_outlets == 0) {
    std::cerr << "Error outlet 1d data should not be empty.\n";
    exit(EXIT_FAILURE);
  }

  auto &input = d_input;

  // step 1: copy relevant data
  for (const auto &a : data_1d) {
    d_perf_pts.push_back(lm::Point(a.p.x, a.p.y, a.p.z));
    d_perf_radii.push_back(a.radius_last);
    d_perf_pres.push_back(a.pressure);
    d_perf_pres_vein.push_back(6666.); // FIXME
    d_perf_nut.push_back(1.); // FIXME
    d_perf_nut_vein.push_back(0.); // FIXME
    d_perf_ball_radii.push_back(0.);
    d_perf_coeff_a.push_back(0.);
    d_perf_coeff_b.push_back(0.);
    d_perf_p_3d_weighted.push_back(0.);
    d_perf_nut_3d_weighted.push_back(0.);
  }

  //
  std::vector<double> perf_flow_capacity;
  for (const auto &r : d_perf_radii)
    perf_flow_capacity.push_back(std::pow(r, 3));
  double max_r3 = max(perf_flow_capacity);
  double min_r3 = min(perf_flow_capacity);

  // step 2: setup perfusion neighborhood
  // instead of point source, we have volume source supported over a ball.
  // radius of ball is proportional to the outlet radius^3 and varies from [ball_r_min, ball_r_max]
  double ball_r_min = input.d_perf_neigh_size.first;  // * d_input.d_h; // avoid point sources
  double ball_r_max = input.d_perf_neigh_size.second; // * d_input.d_h; // avoid too large neighborhood
  std::cout << "ball r = ";
  for (size_t i = 0; i < num_perf_outlets; i++) {
    d_perf_ball_radii[i] = ball_r_min + (ball_r_max - ball_r_min) * (perf_flow_capacity[i] - min_r3) / (max_r3 - min_r3);
    std::cout << d_perf_ball_radii[i] << "; ";
  }
  std::cout << std::endl;

  //  create outlet functions (we choose linear \phi(r) = 1 - r
  for (size_t i = 0; i < num_perf_outlets; i++) {
    if (input.d_perf_fn_type == "const")
      d_perf_fns.push_back(std::make_unique<ConstOutletRadial>(d_perf_pts[i], d_perf_ball_radii[i]));
    else if (input.d_perf_fn_type == "linear")
      d_perf_fns.push_back(std::make_unique<LinearOutletRadial>(d_perf_pts[i], d_perf_ball_radii[i]));
    else if (input.d_perf_fn_type == "gaussian")
      d_perf_fns.push_back(std::make_unique<GaussianOutletRadial>(d_perf_pts[i], d_perf_ball_radii[i], 0.5 * d_perf_ball_radii[i]));
    else {
      std::cerr << "Error: input flag for outlet weight function is invalid.\n";
      exit(EXIT_FAILURE);
    }
  }

  // step 3: for each outlet, create a list of elements and node affected by outlet source and also create coefficients
  d_perf_elems_3D.resize(num_perf_outlets);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &I_out_fn = d_perf_fns[I];
    std::vector<lm::dof_id_type> I_elems;
    for (const auto &elem : d_mesh.active_local_element_ptr_range()) {
      // check if element is inside the ball
      bool any_node_inside_ball = false;
      for (const auto &node : elem->node_index_range()) {
        const lm::Point nodei = elem->node_ref(node);
        auto dx = d_perf_pts[I] - nodei;
        if (dx.norm() < I_out_fn->d_r - 1.e-10) {
          any_node_inside_ball = true;
          break;
        }
      }

      if (any_node_inside_ball)
        add_unique(I_elems, elem->id());
    } // loop over elems

    d_perf_elems_3D[I] = I_elems;
    std::cout << "nelems (I = " << I << ") = " << I_elems.size() << "\n";
  } // loop over outlets

  // compute normalizing coefficient for each outlet weight function
  // note that elements associated to each outlet is now distributed
  // so we first compute local value and then global
  std::vector<double> local_out_normal_const(num_perf_outlets, 0.);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double c = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);
      // init fe
      d_p_cap.init_fe(elem);
      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        c += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]);
      } // quad point loop
    }   // elem loop
    local_out_normal_const[I] = c;
  } // outlet loop

  // we now need to communicate among all processors to compute the global coefficient at processor 0
  std::vector<double> out_normal_c;
  comm_local_to_global(local_out_normal_const, out_normal_c);

  // last thing is to set the normalizing constant of outlet weight function
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    out_fn_I->set_normalize_const(1. / out_normal_c[I]);
    std::cout << "c (I = " << I << ") = " << (*out_fn_I).d_c << "\n";
  }

  // step 4: compute coefficients that we need to exchange with the network system
  std::vector<double> local_out_coeff_a(num_perf_outlets, 0.);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double a = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);

      // init fe
      d_p_cap.init_fe(elem);

      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        a += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * input.d_Lp_art_cap;
      } // quad point loop
    }   // elem loop

    local_out_coeff_a[I] = a;
  } // outlet loop
  std::cout << std::endl;

  // sum distributed coefficients and sync with all processors
  comm_local_to_global(local_out_coeff_a, d_perf_coeff_a);

  // at this point, all processors must have
  // 1. same d_c for outlet weight function
  // 2. same values of coefficient a

  // to verify that all processor have same values of coefficients a and b and normalizing constant
  if (d_input.d_debug_lvl > 0) {
    std::ofstream of;
    of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", d_input.d_out_dir, d_time, get_comm()->rank()));
    of << "x, y, z, r, ball_r, c, a, b, p_3d_weighted\n";
    for (size_t I = 0; I < num_perf_outlets; I++) {
      auto &out_fn_I = d_perf_fns[I];
      of << d_perf_pts[I](0) << ", " << d_perf_pts[I](1) << ", " << d_perf_pts[I](2) << ", "
         << d_perf_radii[I] << ", " << d_perf_ball_radii[I] << ", "
         << (*out_fn_I).d_c << ", " << d_perf_coeff_a[I] << ", "
         << d_perf_coeff_b[I] << ", " << d_perf_p_3d_weighted[I] << ", " << d_perf_nut_3d_weighted[I] << "\n";
    }
    of.close();
  }
}

void HeartToBreast3DSolver::setup_1d3d_partition(const std::vector<VesselTipCurrentCouplingData> &data_1d) {

  auto num_perf_outlets = data_1d.size();
  if (num_perf_outlets == 0) {
    std::cerr << "Error outlet 1d data should not be empty.\n";
    exit(EXIT_FAILURE);
  }

  auto &input = d_input;

  // step 1: copy relevant data
  for (const auto &a : data_1d) {
    d_perf_pts.push_back(lm::Point(a.p.x, a.p.y, a.p.z));
    d_perf_radii.push_back(a.radius_last);
    d_perf_pres.push_back(a.pressure);
    d_perf_pres_vein.push_back(6666.); // FIXME
    d_perf_nut.push_back(1.); // FIXME
    d_perf_nut_vein.push_back(0.); // FIXME
    d_perf_ball_radii.push_back(0.);
    d_perf_coeff_a.push_back(0.);
    d_perf_coeff_b.push_back(0.);
    d_perf_p_3d_weighted.push_back(0.);
    d_perf_nut_3d_weighted.push_back(0.);
  }

  //  create outlet functions (for this case, it will be constant function)
  for (size_t i = 0; i < num_perf_outlets; i++) {
    if (input.d_perf_fn_type == "const")
      d_perf_fns.push_back(std::make_unique<ConstOutletRadial>(d_perf_pts[i], DBL_MAX)); // so that it is practically 1 for any point
    else {
      std::cerr << "Error: input flag for outlet weight function is invalid for partitioned perfusion.\n";
      exit(EXIT_FAILURE);
    }
  }

  // step 3: for each outlet, create a list of elements and node affected by outlet source and also create coefficients
  d_perf_elems_3D.resize(num_perf_outlets);
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {
    auto xc = elem->centroid();

    // find the closest outlet to the center of this element
    double dist = (xc - d_perf_pts[0]).norm();
    long I_found = 0;
    for (size_t I = 0; I < num_perf_outlets; I++) {
      double dist2 = (xc - d_perf_pts[I]).norm();
      if (dist2 < dist) {
        I_found = I;
        dist = dist2;
      }
    }

    // add element to outlet
    d_perf_elems_3D[I_found].push_back(elem->id());
  }

  for (size_t I = 0; I < num_perf_outlets; I++)
    std::cout << "nelems (I = " << I << ") = " << d_perf_elems_3D[I].size() << "\n";

  // compute normalizing coefficient for each outlet weight function
  // note that elements associated to each outlet is now distributed
  // so we first compute local value and then global
  std::vector<double> local_out_normal_const(num_perf_outlets, 0.);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double c = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      c += elem->volume();
    } // elem loop
    local_out_normal_const[I] = c;
  } // outlet loop

  // we now need to communicate among all processors to compute the global coefficient at processor 0
  std::vector<double> out_normal_c;
  comm_local_to_global(local_out_normal_const, out_normal_c);

  // last thing is to set the normalizing constant of outlet weight function
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    out_fn_I->set_normalize_const(1. / out_normal_c[I]);
    std::cout << "c (I = " << I << ") = " << (*out_fn_I).d_c << "\n";
  }

  // step 4: compute coefficients that we need to exchange with the network system
  std::vector<double> local_out_coeff_a(num_perf_outlets, 0.);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double a = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);

      // init fe
      d_p_cap.init_fe(elem);

      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        a += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * input.d_Lp_art_cap;
      } // quad point loop
    }   // elem loop

    local_out_coeff_a[I] = a;
  } // outlet loop

  // sum distributed coefficients and sync with all processors
  comm_local_to_global(local_out_coeff_a, d_perf_coeff_a);

  // at this point, all processors must have
  // 1. same d_c for outlet weight function
  // 2. same values of coefficient a

  // to verify that all processor have same values of coefficients a and b and normalizing constant
  if (d_input.d_debug_lvl > 0) {
    std::ofstream of;
    of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", d_input.d_out_dir, d_time, get_comm()->rank()));
    of << "x, y, z, r, ball_r, c, a, b, p_3d_weighted\n";
    for (size_t I = 0; I < num_perf_outlets; I++) {
      auto &out_fn_I = d_perf_fns[I];
      of << d_perf_pts[I](0) << ", " << d_perf_pts[I](1) << ", " << d_perf_pts[I](2) << ", "
         << d_perf_radii[I] << ", " << d_perf_ball_radii[I] << ", "
         << (*out_fn_I).d_c << ", " << d_perf_coeff_a[I] << ", "
         << d_perf_coeff_b[I] << ", " << d_perf_p_3d_weighted[I] << ", " << d_perf_nut_3d_weighted[I] << "\n";
    }
    of.close();
  }
}

void HeartToBreast3DSolver::comm_local_to_global(const std::vector<double> &local, std::vector<double> &global) {

  auto n = local.size();

  // collect data of all processors in processor 0
  const auto &comm = get_comm();
  std::vector<double> recv_c = local;
  comm->gather(0, recv_c);

  // in rank 0, compute coefficients
  if (comm->rank() == 0) {
    // resize of appropriate size
    global.resize(n);
    for (auto &a : global)
      a = 0.;

    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < n; I++)
        global[I] += recv_c[i * n + I];
    }
  } else
    // in rank other than 0, set the data size to zero so that 'allgather' works properly
    global.resize(0);

  // do allgather (since rank \neq 0 has no data and only rank 0 has data, this should work)
  comm->allgather(global);
}
std::vector<VesselTipCurrentCouplingData3D> HeartToBreast3DSolver::get_vessel_tip_data_3d() {
  std::vector<VesselTipCurrentCouplingData3D> data;
  for (size_t i = 0; i < d_perf_pts.size(); i++) {
    VesselTipCurrentCouplingData3D d;
    d.d_a = d_perf_coeff_a[i];
    d.d_b = d_perf_coeff_b[i];
    d.d_p_3d_w = d_perf_p_3d_weighted[i];
    d.d_nut_3d_w = d_nut_convert_factor * d_perf_nut_3d_weighted[i];
    data.push_back(d);
  }

  return data;
}

void HeartToBreast3DSolver::update_3d_data() {
  auto num_perf_outlets = d_perf_pts.size();

  auto &input = d_input;

  // recompute coefficients b and avg 3d pressure
  std::vector<double> local_out_coeff_b(num_perf_outlets, 0.);
  std::vector<double> local_out_p_3d_weighted(num_perf_outlets, 0.);
  std::vector<double> local_out_nut_3d_weighted(num_perf_outlets, 0.);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double b = 0.;
    double p_3d_w = 0.; // 3D weighted pressure at outlet
    double nut_3d_w = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);
      d_nut_cap.init_dof(elem);

      // init fe
      d_p_cap.init_fe(elem);

      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        // get pressure at quad point
        double p_qp = 0.;
        double nut_qp = 0.;
        for (unsigned int l = 0; l < d_p_cap.d_phi.size(); l++) {
          p_qp += d_p_cap.d_phi[l][qp] * d_p_cap.get_current_sol(l);
          nut_qp += d_p_cap.d_phi[l][qp] * d_nut_cap.get_current_sol(l);
        }

        b += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * input.d_Lp_art_cap * p_qp;

        p_3d_w += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * p_qp;
        nut_3d_w += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * nut_qp;
      } // quad point loop
    }   // elem loop

    local_out_coeff_b[I] = b;
    local_out_p_3d_weighted[I] = p_3d_w;
    local_out_nut_3d_weighted[I] = nut_3d_w;
  } // outlet loop

  // sum distributed coefficients and sync with all processors
  comm_local_to_global(local_out_coeff_b, d_perf_coeff_b);
  comm_local_to_global(local_out_p_3d_weighted, d_perf_p_3d_weighted);
  comm_local_to_global(local_out_nut_3d_weighted, d_perf_nut_3d_weighted);
}

void HeartToBreast3DSolver::update_1d_data(const std::vector<VesselTipCurrentCouplingData> &data_1d) {
  if (d_perf_pts.size() != data_1d.size()) {
    std::cerr << "1D data passed should match the 1D data currently stored in the 3D solver.\n";
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < data_1d.size(); i++) {
    d_perf_pres[i] = data_1d[i].pressure;
    d_perf_pres_vein[i] = 6666.;
    d_perf_nut[i] = data_1d[i].concentration; // FIXME: Check units
    d_perf_nut_vein[i] = 0.;
  }
}

void HeartToBreast3DSolver::set_conductivity_fields() {
  std::vector<unsigned int> dof_indices;
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    d_K_tis_field.get_dof_map().dof_indices(elem, dof_indices);
    d_K_tis_field.solution->set(dof_indices[0], d_input.d_K_tis);

    d_Dnut_tis_field.get_dof_map().dof_indices(elem, dof_indices);
    d_Dnut_tis_field.solution->set(dof_indices[0], d_input.d_Dtis_cap);

    d_N_bar_cap_field.get_dof_map().dof_indices(elem, dof_indices);
    d_N_bar_cap_field.solution->set(dof_indices[0], d_input.d_N_bar_cap);

    d_N_bar_surf_cap_field.get_dof_map().dof_indices(elem, dof_indices);
    d_N_bar_surf_cap_field.solution->set(dof_indices[0], d_input.d_N_bar_surf_cap);
  }

  d_K_tis_field.solution->close();
  d_K_tis_field.update();

  d_Dnut_tis_field.solution->close();
  d_Dnut_tis_field.update();

  d_N_bar_cap_field.solution->close();
  d_N_bar_cap_field.update();

  d_N_bar_surf_cap_field.solution->close();
  d_N_bar_surf_cap_field.update();
}

void HeartToBreast3DSolver::initialize_tumor_field(std::string tumor_mesh_file) {
  d_log("creating tumor mesh\n");
  lm::ReplicatedMesh tum_mesh(*d_comm_p);
  tum_mesh.read(tumor_mesh_file);
  auto tum_mesh_h = get_mesh_size_estimate_using_element_volume(tum_mesh);
  d_log(fmt::format("tumor mesh size = {}\n", tum_mesh_h));

  // create list of tissue mesh element centers
  int nelems = d_mesh.n_elem();
  std::vector<lm::Point> elem_centers(nelems, lm::Point());
  for (const auto &elem : d_mesh.element_ptr_range()) {
    elem_centers[elem->id()] = elem->centroid();
  }

  // create tree for search
  std::unique_ptr<NFlannSearchKd> tree = std::make_unique<NFlannSearchKd>(elem_centers);
  tree->set_input_cloud();

  // loop over elements in tumor mesh and find the closest point in the tree
  for (const auto &tum_elem: tum_mesh.element_ptr_range()) {
    auto xc = tum_elem->centroid();

    // find elements whose center is within radi distance of outlet point
    std::vector<size_t> neighs;
    std::vector<double> sqr_dist;
    auto search_status =
      tree->nearest_search(xc, 5, neighs, sqr_dist);
    for (size_t i=0; i<neighs.size(); i++) {
      if (sqr_dist[i] < 2. * d_input.d_h) {
        // element is initialized with 1 tumor
        size_t tissue_elem_id = neighs[i];
        auto tissue_elem = d_mesh.elem_ptr(tissue_elem_id);

        // set value
        d_tum.init_dof(tissue_elem);
        auto &dof = d_tum.d_dof_indices_sys_var;
        for (size_t j=0; j<dof[0].size(); j++)
          d_tum.d_sys.solution->set(dof[0][j], 1.);
      }
    }
  }
  d_tum.d_sys.solution->close();
  d_tum.d_sys.update();
}

} // namespace macrocirculation

