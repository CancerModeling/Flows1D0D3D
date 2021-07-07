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
#include <memory>

#include "macrocirculation/random_dist.hpp"
#include "macrocirculation/tree_search.hpp"

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"
#include "macrocirculation/mesh_partitioning_perfusion.hpp"
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
  BaseOutletRadial(int type, double r) : d_type(type), d_r(r), d_c(1.) {}

  /* @brief type of radial function. 0 - const, 1 - linear, 2 - quadratic */
  int d_type;

  /* @brief radius of ball on which this function is supported */
  double d_r;

  /* @brief normalizing constant */
  double d_c;

  virtual double operator()(double r) {
    std::cerr << "Error: Function should be defined by inheriting class.\n";
    exit(EXIT_FAILURE);
  }
};

// f(r) = 1, r = [0,1]
class ConstOutletRadial : public BaseOutletRadial {
public:
  ConstOutletRadial(double r) : BaseOutletRadial(0, r) {
    d_c = 3. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(double r) override {
    if (r < d_r - 1.e-10) return d_c;
    else
      return 0.;
  }
};

// f(r) = 1 - r, r = [0,1]
class LinearOutletRadial : public BaseOutletRadial {
public:
  LinearOutletRadial(double r) : BaseOutletRadial(0, r) {
    d_c = 12. / (4. * M_PI * std::pow(d_r, 3)); // normalizing constant
  }

  double operator()(double r) override {
    if (r < d_r - 1.e-10) return d_c * (1. - r / d_r);
    else
      return 0.;
  }
};

// f(r) = exp[-r^2/(2*\sigma^2)], r = [0,1]
class GaussianOutletRadial : public BaseOutletRadial {
public:
  GaussianOutletRadial(double r, double sigma) : BaseOutletRadial(0, r), d_sigma(sigma) {
    double m2 = 0.5 * d_sigma * d_sigma * (d_sigma * std::sqrt(M_PI * 2.) * std::erf(1. / (d_sigma * std::sqrt(2))) - 2. * std::exp(-1. / (2. * d_sigma * d_sigma)));
    d_c = 1. / (4. * M_PI * std::pow(d_r, 3) * m2); // normalizing constant
  }

  /* @brief std of gaussian */
  double d_sigma;

  double operator()(double r) override {
    if (r < d_r - 1.e-10) return d_c * std::exp(-std::pow(r / d_r, 2) / (2. * d_sigma * d_sigma));
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

  void run() override{};

  void write_system(const unsigned int &t_step) override{};

  void solve_system() override{};

  void compute_qoi() override{};

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
  std::vector<std::vector<lm::dof_id_type>> out_nodes;
  std::vector<std::vector<double>> out_coeff;
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
                       lm::ExplicitSystem &hyd_cond,
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

} // namespace darcy3d

int main(int argc, char *argv[]) {
  auto sim_begin = std::chrono::steady_clock::now();

  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value(""))                                               //
    ("final-time", "final simulation time", cxxopts::value<double>()->default_value("1."))                              //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                                    //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.1"))                                          //
    ("hyd-cond", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                               //
    ("permeability", "permeability for mass exchange", cxxopts::value<double>()->default_value("0.1"))                               //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value(""))                                    //
    ("num-points", "number of perfusion points", cxxopts::value<int>()->default_value("10"))                            //
    ("gamma", "value of coefficient for perfusion area estimation", cxxopts::value<double>()->default_value("3"))       //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_darcy3d/")) //
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
  double num_perf_points = args["num-points"].as<int>();

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

  // create outlet point data
  auto &pts = model.pts;
  auto &radii = model.radii;
  auto &out_pres = model.out_pres;

  set_perfusion_pts(out_dir, num_perf_points, pts, radii, hyd_cond, eq_sys, model);
  model.max_r = mc::max(radii);
  model.min_r = mc::min(radii);

  out_pres.resize(num_perf_points);
  mc::DistributionSample<UniformDistribution> uni_dist(10., 1000.);
  for (size_t i = 0; i < num_perf_points; i++)
    out_pres[i] = uni_dist();

  // instead of point source, we have volume source supported over a ball.
  // radius of ball is proportional to the outlet radius and varies from [ball_r_min, ball_r_max]
  double ball_r_min = 0.01 * l;
  double ball_r_max = 0.1 * l;
  auto &ball_r = model.ball_r;
  for (size_t i = 0; i < radii.size(); i++) {
    ball_r.push_back(ball_r_min + (ball_r_max - ball_r_min) * (radii[i] - model.min_r) / (model.max_r - model.min_r));
  }

  // create outlet functions
  auto &out_fns = model.out_fns;
  for (size_t i = 0; i < num_perf_points; i++)
    out_fns.push_back(std::make_unique<darcy3d::LinearOutletRadial>(ball_r[i]));

  // for each outlet, create a list of elements and node affected by outlet source and also create coefficients
  auto &out_elems = model.out_elems;
  auto &out_nodes = model.out_nodes;
  auto &out_coeff = model.out_coeff;

  out_elems.resize(num_perf_points);
  out_nodes.resize(num_perf_points);
  out_coeff.resize(num_perf_points);

  for (size_t I = 0; I < num_perf_points; I++) {

    auto &I_out_fn = out_fns[I];

    std::vector<lm::dof_id_type> I_elems;
    std::vector<lm::dof_id_type> I_nodes;
    std::vector<std::vector<lm::dof_id_type>> I_nodes_elems;

    for (const auto &elem : mesh.active_local_element_ptr_range()) {

      // check if element is inside the ball
      bool any_node_inside_ball = false;
      bool any_node_outside_ball = false;
      for (const auto &node : elem->node_index_range()) {
        const lm::Point nodei = elem->node_ref(node);
        auto dx = pts[I] - nodei;
        if (dx.norm() < ball_r[I] - 1.e-10)
          any_node_inside_ball = true;
        else
          any_node_outside_ball = true;
      }

      if (any_node_inside_ball) {
        mc::add_unique(I_elems, elem->id());
        for (const auto &node : elem->node_index_range()) {
          mc::add_unique(I_nodes, elem->node_ref(node).id());
        }
      }
    } // loop over elems

    // compute coefficients
    auto I_coeff = std::vector<double>(I_nodes.size(), 0.);
    std::vector<unsigned int> dof_indices;
    for (const auto &elem : mesh.active_local_element_ptr_range()) {

      if (mc::locate(I_elems, elem->id()) == -1)
        continue;

      // init dof map
      model.d_p1.init_dof(elem);

      // init fe
      model.d_p1.init_fe(elem);

      for (unsigned int qp = 0; qp < model.d_p1.d_qrule.n_points(); qp++) {

        // compute coefficients
        for (unsigned int i = 0; i < model.d_p1.d_phi.size(); i++) {

          auto node_id = elem->node_id(i);
          auto loc_node = mc::locate(I_nodes, node_id);

          if (loc_node != -1) {
            auto dx = pts[I] - lm::Point(elem->node_ref(i));
            double out_fn_val = (*I_out_fn)(dx.norm());
            I_coeff[loc_node] += model.d_p1.d_JxW[qp] * model.d_p1.d_phi[i][qp] * out_fn_val;
          }
        }
      } // loop over quad points
    }   // loop over elems

    out_nodes[I] = I_nodes;
    out_elems[I] = I_elems;
    out_coeff[I] = I_coeff;
  } // loop over outlets

  // time stepping
  do {
    // Prepare time step
    model.d_step++;
    model.d_time += model.d_dt;

    // change pressure
    mc::DistributionSample<UniformDistribution> uni_dist(10., 1000.);
    for (size_t i = 0; i < num_perf_points; i++)
      out_pres[i] = uni_dist();

    auto solve_clock = std::chrono::steady_clock::now();
    model.d_p1.solve();
    log("pres1 solve time = " + std::to_string(mc::time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

    solve_clock = std::chrono::steady_clock::now();
    model.d_p2.solve();
    log("pres2 solve time = " + std::to_string(mc::time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

    if (model.d_step % 20 == 0) {
      log("writing to file\n");
      mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step / 20) + ".pvtu", eq_sys);
    }
  } while (model.d_time < input.d_T);

  // write
  log("writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step / 20) + ".pvtu", eq_sys);

  return 0;
}

// define assembly functions
void darcy3d::Pres1::assemble() {
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
  auto &out_nodes = d_model_p->out_nodes;
  auto &out_coeff = d_model_p->out_coeff;

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

    // outlet source
    for (size_t I = 0; I < pts.size(); I++) {
      auto &I_out_fn = out_fns[I];
      auto locate_elem = mc::locate(out_elems[I], elem->id());
      if (locate_elem == -1)
        continue;

      for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {

        for (unsigned int i = 0; i < d_phi.size(); i++) {
          if (mc::locate(out_nodes[I], elem->node_id(i)) != -1) {
            double w = (*I_out_fn)((pts[I] - lm::Point(elem->node_ref(i))).norm());
            d_Fe(i) += d_JxW[qp] * d_phi[i][qp] * w;
          }
        }
      } // loop over quad points
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
  auto &out_nodes = d_model_p->out_nodes;
  auto &out_coeff = d_model_p->out_coeff;

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