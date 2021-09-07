////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cxxopts.hpp>
#include <memory>
#include <cstdlib>
#include <fmt/format.h>

#include "macrocirculation/tree_search.hpp"
#include "macrocirculation/random_dist.hpp"

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"
#include "macrocirculation/vtk_io_libmesh.hpp"
#include "macrocirculation/vtk_writer.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

namespace mc = macrocirculation;

// define input, systems, model
namespace darcy3d {
// read and store input parameters
struct InputDeck {
  InputDeck(const std::string &filename = "") : d_K(1.), d_T(1.),
                                                d_dt(0.01), d_h(0.1), d_mesh_file(""),
                                                d_gamma(3), d_nPerfPts(10) {
    if (!filename.empty() and filename != "")
      read_parameters(filename);
  }

  void read_parameters(const std::string &filename) {
    GetPot input(filename);
    d_K = input("K", 1.);
    d_T = input("T", 1.);
    d_dt = input("dt", 0.01);
    d_h = input("h", 0.1);
    d_mesh_file = input("mesh_file", "mesh.e");
    d_gamma = input("gamma", 3);
    d_nPerfPts = input("num_points", 10);
  }

  std::string print_str() {
    std::ostringstream oss;
    oss << "K = " << d_K << "\n";
    oss << "T = " << d_T << "\n";
    oss << "dt = " << d_dt << "\n";
    oss << "h = " << d_h << "\n";
    oss << "mesh_file = " << d_mesh_file << "\n";
    oss << "gamma = " << d_gamma << "\n";
    oss << "num_points = " << d_nPerfPts << "\n";
    return oss.str();
  }

  double d_K;
  double d_T;
  double d_dt;
  double d_h;
  std::string d_mesh_file;
  double d_gamma;
  int d_nPerfPts;
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
class Pres : public mc::BaseAssembly {
public:
  Pres(Model *model, lm::MeshBase &mesh,
       lm::TransientLinearImplicitSystem &sys)
    : mc::BaseAssembly("P", mesh, sys, 1,
                       {sys.variable_number("p")}),
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
        lm::TransientLinearImplicitSystem &pres,
        lm::ExplicitSystem &hyd_cond,
        mc::Logger &log)
    : mc::BaseModel(comm, mesh, eq_sys, log, "Darcy_3D"),
      d_input(input),
      d_pres(this, d_mesh, pres),
      d_hyd_cond(hyd_cond) {
    d_log("model created\n");
  };

  Pres &get_pres() { return d_pres; }

  void run() {};

  void write_system(const unsigned int &t_step) {};

  void solve_system() {};

  void compute_qoi() {};

  InputDeck &d_input;
  Pres d_pres;
  lm::ExplicitSystem &d_hyd_cond;
};

struct NetPoint{
  lm::Point d_x;
  double d_R;
  double d_p;
  double d_Q;
  NetPoint(double gamma = 3.) : d_x(lm::Point()), d_R(0.), d_p(0.), d_Q(std::pow(d_R, gamma)) {}
  NetPoint(lm::Point x, double R, double p = 0., double gamma = 3.) : d_x(x), d_R(R), d_p(p), d_Q(std::pow(d_R, gamma)) {}
};

//
void set_perfusion_pts(std::string out_dir,
                       std::vector<NetPoint> &pts,
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
  int npts = input->d_nPerfPts;
  pts.resize(npts);
  std::vector<int> sel_elems;
  double min_dist = l / npts;
  for (int i=0; i<10*npts; i++) {
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
  mc::DistributionSample<UniformDistribution> uni_dist(min_dist/10., min_dist/3., seed);
  for (int i=0; i<npts; i++)
    pts[i] = NetPoint(elem_centers[sel_elems[i]], uni_dist());
}

//
int create_heterogeneous_conductivity(const lm::MeshBase &mesh,
                                      lm::ExplicitSystem &hyd_cond,
                                      lm::EquationSystems &eq_sys,
                                      const std::vector<int> &elem_taken,
                                      Model &model) {
  std::vector<unsigned int> dof_indices;
  for (const auto &elem : mesh.active_local_element_ptr_range()) {
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);
    int outlet_id = elem_taken[elem->id()];
    hyd_cond.solution->set(dof_indices[0], double(outlet_id+1));
  }
  hyd_cond.solution->close();
  hyd_cond.update();
}

//
void create_perfusion_territory(std::string out_dir,
                                std::vector<NetPoint> &pts,
                                lm::ExplicitSystem &hyd_cond,
                                lm::EquationSystems &eq_sys,
                                Model &model) {

  // initialize random number generator
  int seed = 0;
  srand(seed);

  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  auto &mesh = eq_sys.get_mesh();
  const auto &l = eq_sys.parameters.get<double>("length");

  // create list of element centers for tree search
  int nelems = mesh.n_elem();
  std::vector<lm::Point> elem_centers(nelems, lm::Point());
  double total_vol = 0.;
  for (const auto &elem : mesh.element_ptr_range()) {
    elem_centers[elem->id()] = elem->centroid();
    total_vol += elem->volume();
  }

  // create tree for search
  std::unique_ptr<mc::NFlannSearchKd> tree = std::make_unique<mc::NFlannSearchKd>(elem_centers);
  tree->set_input_cloud();

  // step 1 - assign volume to each outlet point
  model.d_log("  Assigning volumes to each outlet point\n");
  std::vector<double> pts_vol;
  double flow_sum = 0.;
  double min_vol = total_vol;
  for (auto i : pts) flow_sum += i.d_Q;
  for (int i=0; i<pts.size(); i++) {
    double voli = (pts[i].d_Q / flow_sum) * total_vol;
    pts_vol.push_back(voli);

    if (voli < min_vol)
      min_vol = voli;

    model.d_log(fmt::format("    i = {}, voli = {}\n", i, voli));
  }

  double vol_tol = 0.001 * min_vol;

  // step 2 - sort the points based on increasing volume
  model.d_log("  Sorting points\n");
  auto sort_indices = mc::sort_indexes(pts_vol);

  // print
  {
    model.d_log("    info before sorting\n");
    std::ostringstream oss;
    std::vector<size_t> temp;
    for (size_t i=0; i<pts.size(); i++) temp.push_back(i);
    oss << "    pts ids = " << mc::print_str(temp);
    oss << "    vol target = " << mc::print_str(pts_vol);
    model.d_log(oss.str());
  }

  // rearrange
  std::vector<NetPoint> pts_temp = pts;
  std::vector<double> pts_vol_temp = pts_vol;
  std::vector<double> pts_vol_rad(pts.size(), 0.);
  for (size_t i=0; i<sort_indices.size(); i++) {
    auto ii = sort_indices[i];
    pts[i] = pts_temp[ii];
    pts_vol[i] = pts_vol_temp[ii];
    pts_vol_rad[i] = 1.25 * std::pow(3.*pts_vol[i]/(4. * M_PI), 1./3.);
  }

  // print
  {
    model.d_log("    info after sorting\n");
    std::ostringstream oss;
    oss << "    pts ids = " << mc::print_str(sort_indices);
    oss << "    vol target = " << mc::print_str(pts_vol);
    oss << "    radius target = " << mc::print_str(pts_vol_rad);
    model.d_log(oss.str());

    // output vascular domain elements
    auto vtu_writer = mc::VTKWriter(out_dir + "perfusion_points.vtu");

    // create point and radius data
    std::vector<lm::Point> pts_x;
    std::vector<double> pts_r;
    std::vector<double> pts_i;
    for (size_t i=0; i<pts.size(); i++) {
      auto pti = pts[i];
      pts_x.push_back(pti.d_x);
      pts_r.push_back(pti.d_R);
      pts_i.push_back(i+1);
    }
    mc::add_points(pts_x, vtu_writer.d_d_p);
    mc::add_array("Radius", pts_r, vtu_writer.d_d_p);
    mc::add_array("Indices", pts_i, vtu_writer.d_d_p);
    vtu_writer.write();
  }


  // step 3 - brute force element assigning
  model.d_log("  Assigning elements to outlet points using brute-force\n");
  std::vector<std::vector<long>> pts_elem(pts.size(), std::vector<long>());
  std::vector<double> pts_actual_vol(pts.size(), 0.);
  std::vector<int> elem_taken(mesh.n_elem(), -1);
  for (int i=0; i < pts.size(); i++) {
    auto pti = pts[i];
    auto voli = pts_vol[i];
    auto radi = pts_vol_rad[i];

    // find elements whose center is within radi distance of outlet point
    std::vector<size_t> neighs;
    std::vector<double> sqr_dist;
    auto search_status =
      tree->radius_search(pti.d_x, radi, neighs, sqr_dist);

    model.d_log("    search_status = " + std::to_string(search_status) + "\n");

    if (search_status == 0) {
      std::cerr << "Error: Not enough elements in the neighborhood of outlet point\n";
      exit(EXIT_FAILURE);
    }

    double sum_vol_i = 0.;
    for (int j=0; j<neighs.size(); j++) {
      auto j_id = neighs[j];

      // check if this element is already taken
      if (elem_taken[j_id] != -1)
        continue;

      auto elemj_vol = mesh.elem_ptr(j_id)->volume();
      if (sum_vol_i + elemj_vol < voli + vol_tol) {

        pts_elem[i].push_back(j_id);
        sum_vol_i += elemj_vol;

        // store the fact that j_id element now is owned by outlet point i
        elem_taken[j_id] = i;
      }
    }
    pts_actual_vol[i] = sum_vol_i;
  }

  // update conductivity parameter and write to file
  create_heterogeneous_conductivity(mesh, hyd_cond, eq_sys, elem_taken, model);
  model.d_log("  writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(0) + ".pvtu", eq_sys);

  // step 4 - randomly readjust the element distribution
  model.d_log("  Readjustment of element distribution\n");
  std::vector<int> pts_id_vec;
  for (int i=0; i<pts.size(); i++) pts_id_vec.push_back(i);

  int num_random_loop = 200;
  int out_count = 1;
  std::vector<double> store_err;
  std::vector<double> count_no_outlet_elem;
  std::vector<double> vol_not_dist;
  for (int rloop = 0; rloop < num_random_loop; rloop++) {
    model.d_log("    random loop = " + std::to_string(rloop) + "\n");

    // shuffle pts_id_vec
    std::mt19937 gen( std::chrono::system_clock::now().time_since_epoch().count() );
    std::vector<int> v(pts_id_vec.begin(), pts_id_vec.end());
    std::shuffle(v.begin(), v.end(), gen);
    pts_id_vec.assign(v.begin(), v.end());

    std::ostringstream oss;
    oss << "      shuffled list = " << mc::print_str(pts_id_vec);

    std::vector<size_t> temp;
    for (size_t i=0; i<pts.size(); i++) temp.push_back(i);
    oss << "      pts ids = " << mc::print_str(temp);

    temp.clear();
    for (auto i : pts_elem) temp.push_back(i.size());
    oss << "      num elements = " << mc::print_str(temp);

    oss << "      vol actual = " << mc::print_str(pts_actual_vol);
    oss << "      vol target = " << mc::print_str(pts_vol);
    oss << "      radius target = " << mc::print_str(pts_vol_rad);

    std::vector<double> temp2;
    for (int i=0; i<pts.size(); i++) temp2.push_back(100. * std::abs(pts_vol[i] - pts_actual_vol[i])/pts_vol[i]);
    oss << "      percent vol difference = " << mc::print_str(temp2);

    // get average percentage error
    double avg_err = 0.;
    for (auto i : temp2) avg_err += i;
    avg_err = avg_err/temp2.size();
    oss << "      average percent error = " << avg_err << "\n";
    store_err.push_back(avg_err);

    int no_out_el = 0;
    for (auto i : elem_taken) {
      if (i == -1)
        no_out_el++;
    }
    oss << "      number of no outlet elems = " << no_out_el << "\n";
    count_no_outlet_elem.push_back(no_out_el);

    double vol_err = 0.;
    for (auto i : pts_actual_vol) vol_err += i;
    vol_err = 100. * std::abs(total_vol - vol_err)/total_vol;
    oss << "      percent volume not distributed = " << vol_err << "\n";
    vol_not_dist.push_back(vol_err);

    model.d_log(oss.str());


    // loop over outlet points
    int num_pts_change = 0;
    int num_elem_moved = 0;
    for (int ii = 0; ii < pts_id_vec.size(); ii++) {
      int i = pts_id_vec[ii];

      auto pti = pts[i];
      auto voli = pts_vol[i];
      auto radi = pts_vol_rad[i];
      auto voli_act = pts_actual_vol[i];

      std::vector<long> elem_i_old = pts_elem[i];
      auto & elem_i = pts_elem[i];

      // check if we really need to process this outlet
      if (voli_act < voli + vol_tol and voli_act > voli - vol_tol) {
        model.d_log(fmt::format("      skipping outlet end = {}\n", i));
        continue;
      }

      model.d_log(fmt::format("      processing outlet end = {}\n", i));

      // loop over elements of this outlet and for each element find the neighboring element
      int num_neigh_elem = 20;
      double voli_sum_act = voli_act;
      for (auto e : elem_i_old) {
        auto xe = elem_centers[e];

        // if element e is very far from the point center, remove it from the list
        if ((xe - pti.d_x).norm() > 0.5 * l) {
          auto vol_e = mesh.elem_ptr(e)->volume();
          elem_i.erase(std::find(elem_i.begin(), elem_i.end(), e));
          voli_sum_act -= vol_e;

          // free this element
          elem_taken[e] = -1;

          continue;
        }

        std::vector<size_t> ei_neighs;
        std::vector<double> ei_sqr_dist;
        tree->nearest_search(xe, num_neigh_elem, ei_neighs, ei_sqr_dist);

        // loop over neighboring elements
        for (auto ee : ei_neighs) {
          auto vol_ee = mesh.elem_ptr(ee)->volume();
          auto ee_taken = elem_taken[ee];

          // check if ee already exists in the list
          if (ee_taken == i)
            continue;

          // check if adding this element perturbs the volume a lot
          if (voli_sum_act + vol_ee > voli + vol_tol)
            continue;

          // check if at least one adjacent element to this element is owned by outlet i
          bool found_adj_elem_owned_by_i = false;
          auto ee_elem = mesh.elem_ptr(ee);
          for (auto s : ee_elem->side_index_range()) {
            const auto &ee_neigh = ee_elem->neighbor_ptr(s);
            if (ee_neigh != nullptr) {
              found_adj_elem_owned_by_i = elem_taken[ee_neigh->id()] == i;
              if (found_adj_elem_owned_by_i)
                break;
            }
          }

          // if none of the neighbors are owned by i, skip this element
          if (!found_adj_elem_owned_by_i)
            continue;

          // move element ee from previous owner to outlet i
          pts_elem[i].push_back(ee);
          auto ee_outlet_id_old = elem_taken[ee];
          elem_taken[ee] = i;
          voli_sum_act += vol_ee;

          // remove ee from old owner
          if (ee_outlet_id_old != -1) {
            auto &elem_j =pts_elem[ee_outlet_id_old];
            elem_j.erase(std::find(elem_j.begin(), elem_j.end(), ee));
            pts_actual_vol[ee_outlet_id_old] -= vol_ee;
          }
        }
      }

      pts_actual_vol[i] = voli_sum_act;

      if (elem_i_old.size() < pts_elem[i].size()) {
        num_pts_change++;
        num_elem_moved += pts_elem[i].size() -elem_i_old.size();
      }
    }

    model.d_log(fmt::format("      num_pts_change = {}, num_elem_moved = {}\n", num_pts_change, num_elem_moved));

    if (rloop % 10 == 0) {
      // update conductivity parameter and write to file
      create_heterogeneous_conductivity(mesh, hyd_cond, eq_sys, elem_taken, model);
      model.d_log("      writing to file\n");
      mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(out_count) + ".pvtu", eq_sys);
      out_count++;
    }
  }

  {
    std::ofstream oss;
    oss.open(out_dir + "/average_error.txt");
    for (int i=0; i<store_err.size(); i++)
      oss << store_err[i] << ", "
          << count_no_outlet_elem[i] << ", "
          << vol_not_dist[i] << "\n";
    oss.close();
  }
}

} // namespace darcy3d

int main(int argc, char *argv[]) {
  auto sim_begin = std::chrono::steady_clock::now();

  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value(""))                                       //
    ("final-time", "final simulation time", cxxopts::value<double>()->default_value("1."))                      //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                            //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.1"))                                  //
    ("hyd-cond", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                       //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/meshes/darcy_test_tissue_mesh.e"))                            //
    ("gamma", "value of coefficient for perfusion area estimation", cxxopts::value<double>()->default_value("3"))                            //
    ("num-points", "number of perfusion points", cxxopts::value<int>()->default_value("50"))                            //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_part_perfusion_test/")) //
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

  // read input parameters
  auto input = darcy3d::InputDeck(filename);
  if (filename == "") {
    input.d_T = args["final-time"].as<double>();
    input.d_dt = args["time-step"].as<double>();
    input.d_h = args["mesh-size"].as<double>();
    input.d_K = args["hyd-cond"].as<double>();
    input.d_mesh_file = args["mesh-file"].as<std::string>();
    input.d_gamma = args["gamma"].as<double>();
    input.d_nPerfPts = args["num-points"].as<int>();
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
  auto &pres = eq_sys.add_system<lm::TransientLinearImplicitSystem>("P");
  pres.add_variable("p", lm::FIRST);

  // create spatial field of hydraulic conductivity
  auto &hyd_cond = eq_sys.add_system<lm::ExplicitSystem>("K");
  hyd_cond.add_variable("k", lm::CONSTANT, lm::MONOMIAL);

  // create model that holds all essential variables
  log("creating model\n");
  auto model = darcy3d::Model(comm, input, mesh, eq_sys, pres, hyd_cond, log);
  model.d_dt = input.d_dt;

  eq_sys.init();

  // setting up perfusion points
  log("creating random perfusion points\n");
  auto bbox = lm::MeshTools::create_bounding_box(mesh);
  auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
  auto l = (bbox.min() - bbox.max()).norm();
  eq_sys.parameters.set<lm::Point>("center") = xc;
  eq_sys.parameters.set<double>("length") = l;

  std::vector<darcy3d::NetPoint> pts;
  set_perfusion_pts(out_dir, pts, hyd_cond, eq_sys, model);

  // create territory
  create_perfusion_territory(out_dir, pts, hyd_cond, eq_sys, model);

  return 0;
}

// define assembly functions
void darcy3d::Pres::assemble() {}