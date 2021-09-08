////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cxxopts.hpp>
#include <fmt/format.h>
#include <macrocirculation/vtk_writer.hpp>
#include <memory>

#include "macrocirculation/random_dist.hpp"
#include "macrocirculation/tree_search.hpp"

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"
#include "macrocirculation/vtk_io_libmesh.hpp"

namespace mc = macrocirculation;

// define input, systems, model
namespace darcy3d {
// read and store input parameters
struct InputDeck {
  InputDeck(const std::string &filename = "") : d_K(1.), d_T(1.),
                                                d_dt(0.01), d_h(0.1),
                                                d_Lp(0.1), d_mesh_file("") {
    if (!filename.empty() and filename != "")
      read_parameters(filename);
  }

  void read_parameters(const std::string &filename) {
    GetPot input(filename);
    d_K = input("K", 1.);
    d_T = input("T", 1.);
    d_dt = input("dt", 0.01);
    d_h = input("h", 0.1);
    d_Lp = input("Lp", 0.1);
    d_mesh_file = input("mesh_file", "mesh.e");
  }

  std::string print_str() {
    std::ostringstream oss;
    oss << "K = " << d_K << "\n";
    oss << "T = " << d_T << "\n";
    oss << "dt = " << d_dt << "\n";
    oss << "h = " << d_h << "\n";
    oss << "Lp = " << d_Lp << "\n";
    oss << "mesh_file = " << d_mesh_file << "\n";
    return oss.str();
  }

  double d_K;
  double d_T;
  double d_dt;
  double d_h;
  double d_Lp;
  std::string d_mesh_file;
};

class BaseOutletRadial {
public:
  BaseOutletRadial(int type, lm::Point x, double r) : d_type(type), d_x(x), d_r(r), d_c(1.) {}

  /* @brief set the normalizing constant */
  void set_normalize_const(double c) {d_c = c; }

  /* @brief type of radial function. 0 - const, 1 - linear, 2 - quadratic */
  int d_type;

  /* @brief outlet point */
  lm::Point d_x;

  /* @brief radius of ball on which this function is supported */
  double d_r;

  /* @brief normalizing constant */
  double d_c;

  virtual double operator()(const lm::Point &x) const {
    std::cerr << "Error: Function should be defined by inheriting class.\n";
    exit(EXIT_FAILURE);
  }
};

// f(r) = 1, r = [0,1]
class ConstOutletRadial : public BaseOutletRadial {
public:
  ConstOutletRadial(lm::Point x, double r) : BaseOutletRadial(0, x, r) {
    //d_c = 3. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c;
    else
      return 0.;
  }
};

// f(r) = 1 - r, r = [0,1]
class LinearOutletRadial : public BaseOutletRadial {
public:
  LinearOutletRadial(lm::Point x, double r) : BaseOutletRadial(0, x, r) {
    //d_c = 12. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c * (1. - r);
    else
      return 0.;
  }
};

// f(r) = exp[-r^2/(2*\sigma^2)], r = [0,1]
class GaussianOutletRadial : public BaseOutletRadial {
public:
  GaussianOutletRadial(lm::Point x, double r, double sigma) : BaseOutletRadial(0, x, r), d_sigma(sigma) {
    //double m2 = 0.5 * d_sigma * d_sigma * (d_sigma * std::sqrt(M_PI * 2.) * std::erf(1. / (d_sigma * std::sqrt(2))) - 2. * std::exp(-1. / (2. * d_sigma * d_sigma)));
    //d_c = 1. / (4. * M_PI * std::pow(d_r, 3) * m2); // normalizing constant
  }

  /* @brief std of gaussian */
  double d_sigma;

  double operator()(const lm::Point &x) const override {
    double r = (x - d_x).norm() / d_r;
    if (r <= 1. - 1.e-10) return d_c * std::exp(-std::pow(r, 2) / (2. * d_sigma * d_sigma));
    else
      return 0.;
  }
};

// forward declare model class (specific definition requires first defining systems)
class Model;

// define pressure system

// bc
void bc_p(lm::EquationSystems &es) {
  std::set<lm::boundary_id_type> ids;
  ids.insert(2);
  auto &sys = es.get_system<lm::TransientLinearImplicitSystem>("P");
  std::vector<unsigned int> vars(1, sys.variable_number("p"));

  lm::ConstFunction<lm::Number> cf(1.);
  lm::DirichletBoundary bc(ids, vars, &cf);
  sys.get_dof_map().add_dirichlet_boundary(bc);
}

// ic
lm::Number ic_p(const lm::Point &p, const lm::Parameters &es,
                const std::string &system_name,
                const std::string &var_name) { return 0.; }

void ic(lm::EquationSystems &es, const std::string &system_name) {
  auto &sys = es.get_system<lm::TransientLinearImplicitSystem>(
    system_name);
  if (system_name == "P")
    sys.project_solution(ic_p, nullptr,
                         es.parameters);
}

// assembly
class Pres1 : public mc::BaseAssembly {
public:
  Pres1(Model *model, lm::MeshBase &mesh,
        lm::TransientLinearImplicitSystem &sys)
      : mc::BaseAssembly("P1", mesh, sys, 1,
                         {sys.variable_number("p1")}),
        d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
    //bc_p(sys.get_equation_systems()); // add bc
  }

  void assemble() override;
  void assemble_1d();
  Model *d_model_p;
};

class Pres2 : public mc::BaseAssembly {
public:
  Pres2(Model *model, lm::MeshBase &mesh,
        lm::TransientLinearImplicitSystem &sys)
      : mc::BaseAssembly("P2", mesh, sys, 1,
                         {sys.variable_number("p2")}),
        d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
    //bc_p(sys.get_equation_systems()); // add bc
  }

  void assemble() override;
  void assemble_1d();
  Model *d_model_p;
};

// complete definition of model
class Model : public mc::BaseModel {
public:
  Model(lm::Parallel::Communicator *comm,
        InputDeck &input,
        lm::ReplicatedMesh &mesh,
        lm::EquationSystems &eq_sys,
        lm::TransientLinearImplicitSystem &p1,
        lm::TransientLinearImplicitSystem &p2,
        lm::ExplicitSystem &hyd_cond,
        mc::Logger &log)
      : mc::BaseModel(comm, mesh, eq_sys, log, "Darcy_3D"),
        d_input(input),
        d_p1(this, d_mesh, p1),
        d_p2(this, d_mesh, p2),
        d_hyd_cond(hyd_cond) {
    d_log("model created\n");
  };

  void run() {};

  void write_system(const unsigned int &t_step) {};

  void solve_system() {};

  void compute_qoi() {};

  InputDeck &d_input;
  Pres1 d_p1;
  Pres2 d_p2;
  lm::ExplicitSystem &d_hyd_cond;

  std::vector<lm::Point> pts;
  std::vector<double> radii;
  std::vector<double> ball_r;
  std::vector<double> out_pres;
  double max_r;
  double min_r;
  std::vector<std::unique_ptr<darcy3d::BaseOutletRadial>> out_fns;
  std::vector<std::vector<lm::dof_id_type>> out_elems;
  std::vector<double> out_coeff_a;
  std::vector<double> out_coeff_b;
};

//
void create_heterogeneous_conductivity(lm::MeshBase &mesh, lm::ExplicitSystem &hyd_cond, lm::EquationSystems &eq_sys) {

  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  const auto &xc = eq_sys.parameters.get<lm::Point>("center");
  const auto &l = eq_sys.parameters.get<double>("length");

  std::vector<unsigned int> dof_indices;
  for (const auto &elem : mesh.active_local_element_ptr_range()) {
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);
    auto x = elem->centroid();
    //    if ((x - xc).norm() < 0.1 * l)
    //      hyd_cond.solution->set(dof_indices[0], 0.1 * input->d_K);
    //    else
    hyd_cond.solution->set(dof_indices[0], input->d_K);
  }
  hyd_cond.solution->close();
  hyd_cond.update();
}

//
void set_perfusion_pts(std::string out_dir,
                       int num_pts,
                       std::vector<lm::Point> &pts,
                       std::vector<double> &radii,
                       lm::EquationSystems &eq_sys,
                       Model &model) {

  // initialize random number generator
  int seed = 0;
  srand(seed);

  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  const auto &mesh = eq_sys.get_mesh();
  const auto &l = eq_sys.parameters.get<double>("length");

  // create list of element centers for tree search
  int nelems = mesh.n_elem();
  std::vector<lm::Point> elem_centers(nelems, lm::Point());
  for (const auto &elem : mesh.element_ptr_range())
    elem_centers[elem->id()] = elem->centroid();

  // randomly selected desired number of element centers as outlet perfusion points
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
    if (mc::locate_in_set(e, sel_elems) != -1)
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
  mc::DistributionSample<UniformDistribution> uni_dist(min_dist / 10., min_dist / 3., seed);
  for (int i = 0; i < npts; i++) {
    pts[i] = elem_centers[sel_elems[i]];
    radii[i] = uni_dist();
  }
}

void estimate_flow_rate(const std::vector<double> &radii,
                                          double gamma,
                                          std::vector<double> &flow) {
  flow.clear();
  for (size_t i = 0; i < radii.size(); i++)
    flow.push_back(std::pow(radii[i], gamma));
}

void estimate_vol_fraction(const std::vector<double> &radii,
                                             const std::vector<double> &flow,
                                             std::vector<double> &vol) {
  vol.clear();

  double flow_sum = 0.;
  for (auto q : flow) flow_sum += q;
  for (int i = 0; i < flow.size(); i++)
    vol.push_back(flow[i] / flow_sum);
}

// output perfusion points
void output_perfusion_pts(std::string out_file,
                          std::vector<lm::Point> &pts,
                          std::vector<double> &radii,
                          std::vector<double> &ball_r,
                          std::vector<double> &pres,
                          Model &model) {

  // output vascular domain elements
  auto vtu_writer = mc::VTKWriter(out_file);

  mc::add_points(pts, vtu_writer.d_d_p);
  mc::add_array("Radius", radii, vtu_writer.d_d_p);
  mc::add_array("Ball_Radius", ball_r, vtu_writer.d_d_p);
  mc::add_array("pv", pres, vtu_writer.d_d_p);
  vtu_writer.write();
}

void create_outlets(std::string out_dir,
                    int num_perf_points,
                    lm::EquationSystems &eq_sys,
                    Model &model) {

  const auto &mesh = model.get_mesh();
  auto &pres = model.d_p1;

  //
  //  Step 1: setup perfusion outlets and its properties
  //
  // create outlet point data
  auto &pts = model.pts;
  auto &radii = model.radii;
  auto &out_pres = model.out_pres;

  darcy3d::set_perfusion_pts(out_dir, num_perf_points, pts, radii, eq_sys, model);
  model.max_r = mc::max(radii);
  model.min_r = mc::min(radii);

  // random pressure between 10 and 1000
  out_pres.resize(num_perf_points);
  mc::DistributionSample<UniformDistribution> uni_dist(10., 1000., 0);
  for (size_t i = 0; i < num_perf_points; i++)
    out_pres[i] = uni_dist();

  // instead of point source, we have volume source supported over a ball.
  // radius of ball is proportional to the outlet radius and varies from [ball_r_min, ball_r_max]
  double ball_r_min = 4 * model.d_input.d_h;
  double ball_r_max = 10. * model.d_input.d_h;
  auto &ball_r = model.ball_r;
  for (size_t i = 0; i < radii.size(); i++) {
    ball_r.push_back(ball_r_min + (ball_r_max - ball_r_min) * (radii[i] - model.min_r) / (model.max_r - model.min_r));
  }

  //
  //  Step 2: setup perfusion outlet element list and weight function
  //
  // create outlet functions (we choose linear \phi(r) = 1 - r
  auto &out_fns = model.out_fns;
  for (size_t i = 0; i < num_perf_points; i++)
    out_fns.push_back(std::make_unique<darcy3d::LinearOutletRadial>(pts[i], ball_r[i]));

  // for each outlet, create a list of elements and node affected by outlet source and also create coefficients
  auto &out_elems = model.out_elems;
  out_elems.resize(num_perf_points);
  for (size_t I = 0; I < num_perf_points; I++) {
    auto &I_out_fn = out_fns[I];
    std::vector<lm::dof_id_type> I_elems;
    for (const auto &elem : mesh.active_local_element_ptr_range()) {
      // check if element is inside the ball
      bool any_node_inside_ball = false;
      for (const auto &node : elem->node_index_range()) {
        const lm::Point nodei = elem->node_ref(node);
        auto dx = pts[I] - nodei;
        if (dx.norm() < ball_r[I] - 1.e-10) {
          any_node_inside_ball = true;
          break;
        }
      }

      if (any_node_inside_ball)
        mc::add_unique(I_elems, elem->id());
    } // loop over elems

    out_elems[I] = I_elems;
  } // loop over outlets

  // compute coefficients now
  std::vector<double> local_out_normal_const(num_perf_points, 0.);
  for (size_t I=0; I<num_perf_points; I++) {
    auto &out_fn_I = out_fns[I];
    double c = 0.;
    // loop over elements
    for (const auto &elem_id : out_elems[I]) {
      const auto &elem = mesh.elem_ptr(elem_id);
      // init dof map
      pres.init_dof(elem);
      // init fe
      pres.init_fe(elem);
      // loop over quad points
      for (unsigned int qp = 0; qp < pres.d_qrule.n_points(); qp++) {
        c += pres.d_JxW[qp] * (*out_fn_I)(pres.d_qpoints[qp]);
      } // quad point loop
    } // elem loop
    local_out_normal_const[I] = c;
  } // outlet loop

  // we now need to communicate among all processors to compute the total coefficients at processor 0
  const auto &comm = model.get_comm();
  std::vector<double> recv_c = local_out_normal_const;
  comm->gather(0, recv_c);

  // this is temporary for storing normalizing constants
  std::vector<double> out_normal_c;

  // in rank 0, compute coefficients
  if (comm->rank() == 0) {
    // resize of appropriate size
    out_normal_c.resize(num_perf_points);
    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < num_perf_points; I++)
        out_normal_c[I] += recv_c[i*num_perf_points + I];
    }
  }
  else
    // in rank other than 0, get coefficients
    out_normal_c.resize(0);

  // do allgather (since rank/= 0 has no data and only rank =0 has data, this should work)
  comm->allgather(out_normal_c);

  // last thing is to set the normalizing constant
  for (size_t I = 0; I < num_perf_points; I++) {
    auto &out_fn_I = out_fns[I];
    (*out_fn_I).d_c += out_normal_c[I];
  }

  //
  // step 3: compute coefficients that we need to exchange with the network system
  //
  std::vector<double> local_out_coeff_a(num_perf_points, 0.);
  std::vector<double> local_out_coeff_b(num_perf_points, 0.);
  for (size_t I=0; I<num_perf_points; I++) {
    auto &out_fn_I = out_fns[I];
    double a = 0.;
    double b = 0.;
    // loop over elements
    for (const auto &elem_id : out_elems[I]) {
      const auto &elem = mesh.elem_ptr(elem_id);
      // init dof map
      pres.init_dof(elem);
      // init fe
      pres.init_fe(elem);
      // loop over quad points
      for (unsigned int qp = 0; qp < pres.d_qrule.n_points(); qp++) {
        a += pres.d_JxW[qp] * (*out_fn_I)(pres.d_qpoints[qp]) * model.d_input.d_Lp;

        // get pressure at quad point
        double p_qp = 0.;
        for (unsigned int l = 0; l < pres.d_phi.size(); l++) {
          p_qp += pres.d_phi[l][qp] * pres.get_current_sol(l);
        }

        b += pres.d_JxW[qp] * (*out_fn_I)(pres.d_qpoints[qp]) * model.d_input.d_Lp * p_qp;
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
    model.out_coeff_a.resize(num_perf_points);
    model.out_coeff_b.resize(num_perf_points);
    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < num_perf_points; I++) {
        model.out_coeff_a[I] += recv_a[i*num_perf_points + I];
        model.out_coeff_b[I] += recv_b[i*num_perf_points + I];
      }
    }
  }
  else {
    // in rank other than 0, get coefficients
    model.out_coeff_a.resize(0);
    model.out_coeff_b.resize(0);
  }

  // do allgather (since rank/= 0 has no data and only rank =0 has data, this should work)
  comm->allgather(model.out_coeff_a);
  comm->allgather(model.out_coeff_b);

  // at this point, all processors must have
  // 1. same d_c for outlet weight function
  // 2. same values of coefficients a and b

  // to verify, store the values in file
  std::ofstream of;
  of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", out_dir, model.d_time, comm->rank()));
  of << "c, a, b\n";
  for (size_t I = 0; I < num_perf_points; I++) {
    auto &out_fn_I = out_fns[I];
    of << (*out_fn_I).d_c << ", " << model.out_coeff_a[I] << ", " << model.out_coeff_b[I] << "\n";
  }
  of.close();
}

void update_out_coeff_b(std::string out_dir, Model &model) {
  const auto &mesh = model.get_mesh();
  auto &pres = model.d_p1;
  const auto &input = model.d_input;
  int num_perf_points = model.pts.size();
  std::vector<double> local_out_coeff_b(num_perf_points, 0.);
  for (size_t I=0; I<num_perf_points; I++) {
    auto &out_fn_I = model.out_fns[I];
    double b = 0.;
    // loop over elements
    for (const auto &elem_id : model.out_elems[I]) {
      const auto &elem = mesh.elem_ptr(elem_id);
      // init dof map
      pres.init_dof(elem);
      // init fe
      pres.init_fe(elem);
      // loop over quad points
      for (unsigned int qp = 0; qp < pres.d_qrule.n_points(); qp++) {
        // get pressure at quad point
        double p_qp = 0.;
        for (unsigned int l = 0; l < pres.d_phi.size(); l++) {
          p_qp += pres.d_phi[l][qp] * pres.get_current_sol(l);
        }

        b += pres.d_JxW[qp] * (*out_fn_I)(pres.d_qpoints[qp]) * model.d_input.d_Lp * p_qp;
      } // quad point loop
    } // elem loop

    local_out_coeff_b[I] = b;
  } // outlet loop

  // we now need to communicate among all processors to compute the total coefficients at processor 0
  const auto &comm = model.get_comm();
  std::vector<double> recv_b = local_out_coeff_b;
  comm->gather(0, recv_b);
  // in rank 0, compute coefficients
  if (comm->rank() == 0) {
    // resize of appropriate size
    model.out_coeff_b.resize(num_perf_points);
    for (size_t i=0; i<num_perf_points; i++)
      model.out_coeff_b[i] = 0.;
    // compute
    for (int i = 0; i < comm->size(); i++) {
      for (size_t I = 0; I < num_perf_points; I++) {
        model.out_coeff_b[I] += recv_b[i*num_perf_points + I];
      }
    }
  }
  else {
    // in rank other than 0, get coefficients
    model.out_coeff_b.resize(0);
  }

  // do allgather (since rank/= 0 has no data and only rank =0 has data, this should work)
  comm->allgather(model.out_coeff_b);

  // to verify, store the values in file
  std::ofstream of;
  of.open(fmt::format("{}outlet_coefficients_t_{:5.3f}_proc_{}.txt", out_dir, model.d_time, comm->rank()));
  of << "c, a, b\n";
  for (size_t I = 0; I < num_perf_points; I++) {
    auto &out_fn_I = model.out_fns[I];
    of << (*out_fn_I).d_c << ", " << model.out_coeff_a[I] << ", " << model.out_coeff_b[I] << "\n";
  }
  of.close();
}

} // namespace darcy3d

int main(int argc, char *argv[]) {
  auto sim_begin = std::chrono::steady_clock::now();

  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value(""))                                                              //
    ("final-time", "final simulation time", cxxopts::value<double>()->default_value("5."))                                             //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                                                   //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.02"))                                                         //
    ("hyd-cond", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                                              //
    ("permeability", "permeability for mass exchange", cxxopts::value<double>()->default_value("0.1"))                                 //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value(""))                                                   //
    ("num-points", "number of perfusion points", cxxopts::value<int>()->default_value("10"))                                           //
    ("gamma", "value of coefficient for perfusion area estimation", cxxopts::value<double>()->default_value("3"))                      //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_darcy_flow_reg_outlets/")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if (!args.unmatched().empty()) {
    std::cout << "The following arguments were unmatched: " << std::endl;
    for (auto &it : args.unmatched())
      std::cout << " " << it;
    std::cout
      << "\nAre they part of petsc or a different auxiliary library?"
      << std::endl;
  }

  auto filename = args["input-file"].as<std::string>();
  auto out_dir = args["output-directory"].as<std::string>();
  double gamma = args["gamma"].as<double>();
  int num_perf_points = args["num-points"].as<int>();

  // read input parameters
  auto input = darcy3d::InputDeck(filename);
  if (filename == "") {
    input.d_T = args["final-time"].as<double>();
    input.d_dt = args["time-step"].as<double>();
    input.d_h = args["mesh-size"].as<double>();
    input.d_K = args["hyd-cond"].as<double>();
    input.d_Lp = args["permeability"].as<double>();
    input.d_mesh_file = args["mesh-file"].as<std::string>();
  }

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  log("input data \n" + input.print_str() + "\n");

  // create mesh
  log("creating mesh\n");
  lm::ReplicatedMesh mesh(*comm);
  long N = long(1. / input.d_h);
  if (input.d_mesh_file != "")
    mesh.read(input.d_mesh_file);
  else
    lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                          1., 0., 1., lm::HEX8);

  // create equation system
  log("creating equation system\n");
  lm::EquationSystems eq_sys(mesh);
  eq_sys.parameters.set<darcy3d::InputDeck *>("input_deck") = &input;
  eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
  auto &p1 = eq_sys.add_system<lm::TransientLinearImplicitSystem>("P1");
  p1.add_variable("p1", lm::FIRST);
  auto &p2 = eq_sys.add_system<lm::TransientLinearImplicitSystem>("P2");
  p2.add_variable("p2", lm::FIRST);

  // create spatial field of hydraulic conductivity
  auto &hyd_cond = eq_sys.add_system<lm::ExplicitSystem>("K");
  hyd_cond.add_variable("k", lm::CONSTANT, lm::MONOMIAL);

  // create model that holds all essential variables
  log("creating model\n");
  auto model = darcy3d::Model(comm, input, mesh, eq_sys, p1, p2, hyd_cond, log);
  model.d_dt = input.d_dt;

  eq_sys.init();

  // create heterogeneous property field
  log("setting up K field\n");
  auto bbox = lm::MeshTools::create_bounding_box(mesh);
  auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
  auto l = (bbox.min() - bbox.max()).norm();
  eq_sys.parameters.set<lm::Point>("center") = xc;
  eq_sys.parameters.set<double>("length") = l;
  darcy3d::create_heterogeneous_conductivity(mesh, hyd_cond, eq_sys);

  // create outlets
  create_outlets(out_dir, num_perf_points, eq_sys, model);

  // time stepping
  do {
    // Prepare time step
    model.d_step++;
    model.d_time += model.d_dt;

    // change pressure
    mc::DistributionSample<UniformDistribution> uni_dist(10., 1000., model.d_step);
    for (size_t i = 0; i < num_perf_points; i++)
      model.out_pres[i] = uni_dist();

    auto solve_clock = std::chrono::steady_clock::now();
    model.d_p1.solve();
    log("pres1 solve time = " + std::to_string(mc::time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

    solve_clock = std::chrono::steady_clock::now();
    //model.d_p2.solve();
    //log("pres2 solve time = " + std::to_string(mc::time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

    if (model.d_step % 20 == 0) {
      log("writing to file\n");
      mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step / 20) + ".pvtu", eq_sys);

      darcy3d::output_perfusion_pts(out_dir + "perfusion_output_" + std::to_string(model.d_step / 20) + ".vtu",
                                    model.pts,
                                    model.radii,
                                    model.ball_r,
                                    model.out_pres,
                                    model);

      // also update the weights
      darcy3d::update_out_coeff_b(out_dir, model);
    }
  } while (model.d_time < input.d_T);

  // write
  log("writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step / 20) + ".pvtu", eq_sys);

  return 0;
}

// define assembly functions
void darcy3d::Pres1::assemble() {

  assemble_1d();

  auto &eq_sys = d_model_p->get_system();
  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  const auto &xc = eq_sys.parameters.get<lm::Point>("center");
  const auto &l = eq_sys.parameters.get<double>("length");
  const double dt = d_model_p->d_dt;
  const double t = d_model_p->d_time;
  auto &hyd_cond = d_model_p->d_hyd_cond;
  std::vector<unsigned int> dof_indices;
  auto &out_fns = d_model_p->out_fns;
  auto &pts = d_model_p->pts;
  auto &ball_r = d_model_p->ball_r;
  auto &out_pres = d_model_p->out_pres;
  auto &out_elems = d_model_p->out_elems;

  lm::Point source_xc = xc + lm::Point(0.1 * l, 0., 0.2 * l);

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);

    // init fe
    init_fe(elem);

    // get K at this element
    double elem_K = hyd_cond.current_solution(dof_indices[0]);

    double rhs = 0.;
    auto x = elem->centroid();
    //if ((x - source_xc).norm() < 0.15 * l)
    //  rhs = std::sin(t / 0.1) * std::exp(-t / 10.);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
      double lhs = d_JxW[qp] * elem_K;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += d_JxW[qp] * rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++)
          d_Ke(i, j) += lhs * d_dphi[j][qp] * d_dphi[i][qp];
      }
    } // loop over quad points

    // dirichlet bc
    {
      // The penalty value.
      const double penalty = 1.e10;

      for (auto s : elem->side_index_range())
        if (elem->neighbor_ptr(s) == nullptr) {
          init_face_fe(elem, s);

          for (unsigned int qp = 0; qp < d_qrule_face.n_points(); qp++) {
            // Matrix contribution
            for (unsigned int i = 0; i < d_phi_face.size(); i++) {
              d_Fe(i) += penalty * d_JxW_face[qp] * d_phi_face[i][qp];
              for (unsigned int j = 0; j < d_phi_face.size(); j++)
                d_Ke(i, j) += penalty * d_JxW_face[qp] * d_phi_face[i][qp] * d_phi_face[j][qp];
            }
          }
        }
    }

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}

void darcy3d::Pres1::assemble_1d() {
  auto &eq_sys = d_model_p->get_system();
  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  auto &out_fns = d_model_p->out_fns;
  auto &pts = d_model_p->pts;
  auto &ball_r = d_model_p->ball_r;
  auto &out_pres = d_model_p->out_pres;
  auto &out_elems = d_model_p->out_elems;

  for (size_t I = 0; I < pts.size(); I++) {
    auto &out_fn_I = out_fns[I];
    const auto &pI = out_pres[I];
    const auto &RI = ball_r[I];
    const auto &xI = pts[I];

    // loop over elements
    for (const auto &elem_id : out_elems[I]) {

      const auto &elem = d_mesh.elem_ptr(elem_id);

      // init dof map
      init_dof(elem);

      // init fe
      init_fe(elem);

      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
        for (unsigned int i = 0; i < d_phi.size(); i++) {
          d_Fe(i) += d_JxW[qp] * input->d_Lp * pI * (*out_fn_I)(d_qpoints[qp]) *d_phi[i][qp];

          for (unsigned int j = 0; j < d_phi.size(); j++)
            d_Ke(i, j) += d_JxW[qp] * input->d_Lp * (*out_fn_I)(d_qpoints[qp]) *d_phi[i][qp] * d_phi[j][qp];
        }
      } // loop over quad points

      d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                       d_dof_indices_sys);
      d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
      d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
    } // elem loop
  }   // outlet loop
}

// define assembly functions
void darcy3d::Pres2::assemble() {
  auto &eq_sys = d_model_p->get_system();
  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  const auto &xc = eq_sys.parameters.get<lm::Point>("center");
  const auto &l = eq_sys.parameters.get<double>("length");
  const double dt = d_model_p->d_dt;
  const double t = d_model_p->d_time;
  auto &hyd_cond = d_model_p->d_hyd_cond;
  double Lp = input->d_Lp;
  std::vector<unsigned int> dof_indices;
  auto &out_fns = d_model_p->out_fns;
  auto &pts = d_model_p->pts;
  auto &ball_r = d_model_p->ball_r;
  auto &out_pres = d_model_p->out_pres;
  auto &out_elems = d_model_p->out_elems;

  lm::Point source_xc = xc + lm::Point(0.1 * l, 0., 0.2 * l);

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);

    // init fe
    init_fe(elem);

    // get K at this element
    double elem_K = hyd_cond.current_solution(dof_indices[0]);

    double rhs = 0.;
    auto x = elem->centroid();
    if ((x - source_xc).norm() < 0.15 * l)
      rhs = std::sin(t / 0.1) * std::exp(-t / 10.);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
      double lhs = d_JxW[qp] * elem_K;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += d_JxW[qp] * rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++)
          d_Ke(i, j) += lhs * d_dphi[j][qp] * d_dphi[i][qp];
      }
    } // loop over quad points

    // dirichlet bc
    {
      // The penalty value.
      const double penalty = 1.e10;

      for (auto s : elem->side_index_range())
        if (elem->neighbor_ptr(s) == nullptr) {
          init_face_fe(elem, s);

          for (unsigned int qp = 0; qp < d_qrule_face.n_points(); qp++) {
            // Matrix contribution
            for (unsigned int i = 0; i < d_phi_face.size(); i++) {
              d_Fe(i) += penalty * d_JxW_face[qp] * d_phi_face[i][qp];
              for (unsigned int j = 0; j < d_phi_face.size(); j++)
                d_Ke(i, j) += penalty * d_JxW_face[qp] * d_phi_face[i][qp] * d_phi_face[j][qp];
            }
          }
        }
    }

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}