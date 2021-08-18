////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "heart_to_breast_3d_solver.hpp"
#include "heart_to_breast_1d_solver.hpp"
#include "vtk_writer.hpp"
#include "random_dist.hpp"
#include "tree_search.hpp"
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
} // local namespace

// input class definitions

HeartToBreast3DSolverInputDeck::HeartToBreast3DSolverInputDeck(const std::string &filename)
    : d_K_cap(1.), d_K_tis(1.), d_Lp_cap(0.1), d_Lp_tis(0.1),
      d_T(1.), d_dt(0.01), d_h(0.1), d_mesh_file(""), d_out_dir(""),
      d_perf_fn_type("linear"), d_perf_neigh_size({4., 10.}),
      d_debug_lvl(0) {
  if (!filename.empty() and filename != "")
    read_parameters(filename);
}

void HeartToBreast3DSolverInputDeck::read_parameters(const std::string &filename) {
  GetPot input(filename);
  d_K_cap = input("K_cap", 1.);
  d_K_tis = input("K_tis", 0.1);
  d_Lp_cap = input("Lp_cap", 1.);
  d_Lp_tis = input("Lp_tis", 0.1);
  d_T = input("T", 1.);
  d_dt = input("dt", 0.01);
  d_h = input("h", 0.1);
  d_mesh_file = input("mesh_file", "");
  d_out_dir = input("out_dir", "");
}

std::string HeartToBreast3DSolverInputDeck::print_str() {
  std::ostringstream oss;
  oss << "K_cap = " << d_K_cap << "\n";
  oss << "K_tis = " << d_K_tis << "\n";
  oss << "Lp_cap = " << d_Lp_cap << "\n";
  oss << "Lp_tis = " << d_Lp_tis << "\n";
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
                                             lm::ExplicitSystem &K_cap_field,
                                             lm::ExplicitSystem &Lp_cap_field,
                                             Logger &log)
    : BaseModel(libmesh_comm, mesh, eq_sys, log, "HeartToBreast3DSolver"),
      d_input(input),
      d_p_cap(this, d_mesh, p_cap),
      d_p_tis(this, d_mesh, p_tis),
      d_K_cap_field(K_cap_field),
      d_Lp_cap_field(Lp_cap_field) {

  d_dt = input.d_dt;
  d_log("created HeartToBreast3DSolver object\n");
}

void HeartToBreast3DSolver::write_perfusion_output(std::string out_file) {

  auto vtu_writer = VTKWriter(out_file);
  add_points(d_perf_pts, vtu_writer.d_d_p);
  add_array("Radius", d_perf_radii, vtu_writer.d_d_p);
  add_array("Ball_Radius", d_perf_ball_radii, vtu_writer.d_d_p);
  add_array("pv", d_perf_pres, vtu_writer.d_d_p);
  vtu_writer.write();
}

void HeartToBreast3DSolver::setup_random_outlets(unsigned int num_perf_outlets) {

  const auto &mesh = get_mesh();

  //  Step 1: setup perfusion outlets and its properties
  set_perfusion_pts(d_input.d_out_dir, num_perf_outlets, d_perf_pts, d_perf_radii, d_eq_sys, *this);
  double max_r = max(d_perf_radii);
  double min_r = min(d_perf_radii);

  // random pressure between 10 and 1000
  d_perf_pres.resize(num_perf_outlets);
  DistributionSample<UniformDistribution> uni_dist(10., 1000., 0);
  for (size_t i = 0; i < num_perf_outlets; i++)
    d_perf_pres[i] = uni_dist();

  // instead of point source, we have volume source supported over a ball.
  // radius of ball is proportional to the outlet radius and varies from [ball_r_min, ball_r_max]
  double ball_r_min = 4 * d_input.d_h;
  double ball_r_max = 10. * d_input.d_h;
  d_perf_pres.clear();
  for (size_t i = 0; i < num_perf_outlets; i++)
    d_perf_pres.push_back(ball_r_min + (ball_r_max - ball_r_min) * (d_perf_radii[i] - min_r) / (max_r - min_r));

  //  Step 2: setup perfusion outlet element list and weight function
  // create outlet functions (we choose linear \phi(r) = 1 - r
  for (size_t i = 0; i < num_perf_outlets; i++)
    d_perf_fns.push_back(std::make_unique<LinearOutletRadial>(d_perf_pts[i], d_perf_pres[i]));

  // for each outlet, create a list of elements and node affected by outlet source and also create coefficients
  d_perf_elems_3D.resize(num_perf_outlets);
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &I_out_fn = d_perf_fns[I];
    std::vector<lm::dof_id_type> I_elems;
    for (const auto &elem : mesh.active_local_element_ptr_range()) {
      // check if element is inside the ball
      bool any_node_inside_ball = false;
      for (const auto &node : elem->node_index_range()) {
        const lm::Point nodei = elem->node_ref(node);
        auto dx = d_perf_pts[I] - nodei;
        if (dx.norm() < d_perf_pres[I] - 1.e-10) {
          any_node_inside_ball = true;
          break;
        }
      }

      if (any_node_inside_ball)
        add_unique(I_elems, elem->id());
    } // loop over elems

    d_perf_elems_3D[I] = I_elems;
  } // loop over outlets

  // compute normalizing coefficient for each outlet weight function
  std::vector<double> local_out_normal_const(num_perf_outlets, 0.);
  for (size_t I=0; I<num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double c = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);
      // init fe
      d_p_cap.init_fe(elem);
      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        c += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]);
      } // quad point loop
    } // elem loop
    local_out_normal_const[I] = c;
  } // outlet loop

  // we now need to communicate among all processors to compute the total coefficient at processor 0
  const auto &comm = get_comm();
  std::vector<double> recv_c = local_out_normal_const;
  comm->gather(0, recv_c);

  // this is temporary for storing normalizing constants
  std::vector<double> out_normal_c;

  // in rank 0, compute coefficients
  if (comm->rank() == 0) {
    // resize of appropriate size
    out_normal_c.resize(num_perf_outlets);
    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < num_perf_outlets; I++)
        out_normal_c[I] += recv_c[i*num_perf_outlets + I];
    }
  }
  else
    // in rank other than 0, get coefficients
    out_normal_c.resize(0);

  // do allgather (since rank/= 0 has no data and only rank =0 has data, this should work)
  comm->allgather(out_normal_c);

  // last thing is to set the normalizing constant
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    (*out_fn_I).d_c += out_normal_c[I];
  }

  // step 3: compute coefficients that we need to exchange with the network system
  std::vector<unsigned int> Lp_cap_dof_indices;
  std::vector<double> local_out_coeff_a(num_perf_outlets, 0.);
  std::vector<double> local_out_coeff_b(num_perf_outlets, 0.);
  for (size_t I=0; I<num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double a = 0.;
    double b = 0.;
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);
      d_Lp_cap_field.get_dof_map().dof_indices(elem, Lp_cap_dof_indices);

      // init fe
      d_p_cap.init_fe(elem);

      // get Lp at this element
      double Lp_elem = d_Lp_cap_field.current_solution(Lp_cap_dof_indices[0]);

      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        a += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * Lp_elem;

        // get pressure at quad point
        double p_qp = 0.;
        for (unsigned int l = 0; l < d_p_cap.d_phi.size(); l++) {
          p_qp += d_p_cap.d_phi[l][qp] * d_p_cap.get_current_sol(l);
        }

        b += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * Lp_elem * p_qp;
      } // quad point loop
    } // elem loop

    local_out_coeff_a[I] = a;
    local_out_coeff_b[I] = b;
  } // outlet loop

  std::vector<double> recv_a = local_out_coeff_a;
  std::vector<double> recv_b = local_out_coeff_b;
  comm->gather(0, recv_a);
  comm->gather(0, recv_b);

  // in rank 0, compute coefficients
  if (comm->rank() == 0) {
    // resize of appropriate size
    d_perf_coeff_a.resize(num_perf_outlets);
    d_perf_coeff_b.resize(num_perf_outlets);
    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < num_perf_outlets; I++) {
        d_perf_coeff_a[I] += recv_a[i*num_perf_outlets + I];
        d_perf_coeff_b[I] += recv_b[i*num_perf_outlets + I];
      }
    }
  }
  else {
    // in rank other than 0, get coefficients
    d_perf_coeff_a.resize(0);
    d_perf_coeff_b.resize(0);
  }

  // do allgather (since rank/= 0 has no data and only rank =0 has data, this should work)
  comm->allgather(d_perf_coeff_a);
  comm->allgather(d_perf_coeff_b);

  // at this point, all processors must have
  // 1. same d_c for outlet weight function
  // 2. same values of coefficients a and b

  // to verify, store the values in file
  std::ofstream of;
  of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", d_input.d_out_dir, d_time, comm->rank()));
  of << "c, a, b\n";
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    of << (*out_fn_I).d_c << ", " << d_perf_coeff_a[I] << ", " << d_perf_coeff_b[I] << "\n";
  }
  of.close();
}

void HeartToBreast3DSolver::setup() {
  // TODO setup other things
  d_eq_sys.init();
}
double HeartToBreast3DSolver::get_time() const {
  return d_time;
}
void HeartToBreast3DSolver::solve() {
}
void HeartToBreast3DSolver::write_output() {
}
void HeartToBreast3DSolver::set_output_folder(std::string output_dir) {
}
void HeartToBreast3DSolver::setup_1d3d(const std::vector<VesselTipCurrentCouplingData> &data_1d) {

  auto num_perf_outlets = data_1d.size();

  // step 1: copy relevant data
  for (const auto &a : data_1d) {
    d_perf_pts.push_back(lm::Point(a.p.x, a.p.y, a.p.z));
    d_perf_radii.push_back(a.radius);
    d_perf_pres.push_back(a.pressure);
    d_perf_ball_radii.push_back(0.);
    d_perf_coeff_a.push_back(0.);
    d_perf_coeff_b.push_back(0.);
    d_perf_p_3d_weighted.push_back(0.);
  }

  //
  std::vector<double> perf_flow_capacity;
  for (const auto & r: d_perf_radii)
    perf_flow_capacity.push_back(std::pow(r, 3));
  double max_r3 = max(perf_flow_capacity);
  double min_r3 = min(perf_flow_capacity);

  // step 2: setup perfusion neighborhood
  // instead of point source, we have volume source supported over a ball.
  // radius of ball is proportional to the outlet radius^3 and varies from [ball_r_min, ball_r_max]
  double ball_r_min = d_input.d_perf_neigh_size.first * d_input.d_h; // avoid point sources
  double ball_r_max = d_input.d_perf_neigh_size.second * d_input.d_h; // avoid too large neighborhood
  for (size_t i = 0; i < num_perf_outlets; i++)
    d_perf_ball_radii[i] = ball_r_min + (ball_r_max - ball_r_min) * (perf_flow_capacity[i] - min_r3) / (max_r3 - min_r3);

  //  create outlet functions (we choose linear \phi(r) = 1 - r
  for (size_t i = 0; i < num_perf_outlets; i++) {
    if (d_input.d_perf_fn_type == "const")
      d_perf_fns.push_back(std::make_unique<ConstOutletRadial>(d_perf_pts[i], d_perf_ball_radii[i]));
    else if (d_input.d_perf_fn_type == "linear")
      d_perf_fns.push_back(std::make_unique<LinearOutletRadial>(d_perf_pts[i], d_perf_ball_radii[i]));
    else if (d_input.d_perf_fn_type == "gaussian")
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
        if (dx.norm() < d_perf_pres[I] - 1.e-10) {
          any_node_inside_ball = true;
          break;
        }
      }

      if (any_node_inside_ball)
        add_unique(I_elems, elem->id());
    } // loop over elems

    d_perf_elems_3D[I] = I_elems;
  } // loop over outlets

  // compute normalizing coefficient for each outlet weight function
  // note that elements associated to each outlet is now distributed
  // so we first compute local value and then global
  std::vector<double> local_out_normal_const(num_perf_outlets, 0.);
  for (size_t I=0; I<num_perf_outlets; I++) {
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
    } // elem loop
    local_out_normal_const[I] = c;
  } // outlet loop

  // we now need to communicate among all processors to compute the global coefficient at processor 0
  std::vector<double> out_normal_c;
  comm_local_to_global(local_out_normal_const, out_normal_c);

  // last thing is to set the normalizing constant of outlet weight function
  for (size_t I = 0; I < num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    (*out_fn_I).d_c += out_normal_c[I];
  }

  // step 4: compute coefficients that we need to exchange with the network system
  std::vector<unsigned int> Lp_cap_dof_indices;
  std::vector<double> local_out_coeff_a(num_perf_outlets, 0.);
  std::vector<double> local_out_coeff_b(num_perf_outlets, 0.);
  std::vector<double> local_out_p_3d_weighted(num_perf_outlets, 0.);
  for (size_t I=0; I<num_perf_outlets; I++) {
    auto &out_fn_I = d_perf_fns[I];
    double a = 0.;
    double b = 0.;
    double p_3d_w = 0.; // 3D weighted pressure at outlet
    // loop over elements
    for (const auto &elem_id : d_perf_elems_3D[I]) {
      const auto &elem = d_mesh.elem_ptr(elem_id);
      // init dof map
      d_p_cap.init_dof(elem);
      d_Lp_cap_field.get_dof_map().dof_indices(elem, Lp_cap_dof_indices);

      // init fe
      d_p_cap.init_fe(elem);

      // get Lp at this element
      double Lp_elem = d_Lp_cap_field.current_solution(Lp_cap_dof_indices[0]);

      // loop over quad points
      for (unsigned int qp = 0; qp < d_p_cap.d_qrule.n_points(); qp++) {
        a += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * Lp_elem;

        // get pressure at quad point
        double p_qp = 0.;
        for (unsigned int l = 0; l < d_p_cap.d_phi.size(); l++) {
          p_qp += d_p_cap.d_phi[l][qp] * d_p_cap.get_current_sol(l);
        }

        b += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * Lp_elem * p_qp;

        p_3d_w += d_p_cap.d_JxW[qp] * (*out_fn_I)(d_p_cap.d_qpoints[qp]) * p_qp;
      } // quad point loop
    } // elem loop

    local_out_coeff_a[I] = a;
    local_out_coeff_b[I] = b;
    local_out_p_3d_weighted[I] = p_3d_w;
  } // outlet loop

  // sum distributed coefficients and sync with all processors
  comm_local_to_global(local_out_coeff_a, d_perf_coeff_a);
  comm_local_to_global(local_out_coeff_b, d_perf_coeff_b);
  comm_local_to_global(local_out_p_3d_weighted, d_perf_p_3d_weighted);

  // at this point, all processors must have
  // 1. same d_c for outlet weight function
  // 2. same values of coefficients a and b

  // to verify that all processor have same values of coefficients a and b and normalizing constant
  if (d_input.d_debug_lvl > 0) {
    std::ofstream of;
    of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", d_input.d_out_dir, d_time, get_comm()->rank()));
    of << "c, a, b, p_3d_weighted\n";
    for (size_t I = 0; I < num_perf_outlets; I++) {
      auto &out_fn_I = d_perf_fns[I];
      of << (*out_fn_I).d_c << ", "
         << d_perf_coeff_a[I] << ", "
         << d_perf_coeff_b[I] << ", "
         << d_perf_p_3d_weighted[I] << "\n";
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
    for (auto & a: global)
      a = 0.;

    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < n; I++)
        global[I] += recv_c[i*n + I];
    }
  }
  else
    // in rank other than 0, set the data size to zero so that 'allgather' works properly
    global.resize(0);

  // do allgather (since rank \neq 0 has no data and only rank 0 has data, this should work)
  comm->allgather(global);
}
std::vector<VesselTipCurrentCouplingData3D> HeartToBreast3DSolver::get_vessel_tip_data_3d() {
  std::vector<VesselTipCurrentCouplingData3D> data;
  for (size_t i=0; i<d_perf_pts.size(); i++) {
    VesselTipCurrentCouplingData3D d;
    d.d_a = d_perf_coeff_a[i];
    d.d_b = d_perf_coeff_b[i];
    d.d_p_3d_w = d_perf_p_3d_weighted[i];
  }

  return data;
}
void HeartToBreast3DSolver::update_1d_data(const std::vector<VesselTipCurrentCouplingData> &data_1d) {
  if (d_perf_pts.size() != data_1d.size()) {
    std::cerr << "1D data passed should match the 1D data currently stored in the 3D solver.\n";
    exit(EXIT_FAILURE);
  }

  for (size_t i=0; i<data_1d.size(); i++)
    d_perf_pres[i] = data_1d[i].pressure;
}

} // namespace macrocirculation