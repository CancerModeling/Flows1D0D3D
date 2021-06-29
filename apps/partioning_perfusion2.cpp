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
#include "macrocirculation/mesh_partitioning_perfusion.hpp"
#include "macrocirculation/vtk_io.hpp"


namespace mc = macrocirculation;

// define input, systems, model
namespace darcy3d {
//
void set_perfusion_pts(std::string out_dir,
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
  int npts = input->d_nPerfPts;
  pts.resize(npts);
  radii.resize(npts);
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
  for (int i=0; i<npts; i++) {
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
                        cxxopts::value<std::string>()->default_value(""))                                       //
    ("final-time", "final simulation time", cxxopts::value<double>()->default_value("1."))                      //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                            //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.1"))                                  //
    ("hyd-cond", "hydraulic conductivity", cxxopts::value<double>()->default_value("1."))                       //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/meshes/darcy_test_tissue_mesh.e"))                            //
    ("gamma", "value of coefficient for perfusion area estimation", cxxopts::value<double>()->default_value("3"))                            //
    ("num-points", "number of perfusion points", cxxopts::value<int>()->default_value("10"))                            //
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
  double gamma = args["gamma"].as<double>();
  double nPerfPts = args["num-points"].as<int>();
  auto mesh_file = args["mesh-file"].as<std::string>();
  auto mesh_size = args["mesh-size"].as<double>();

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  log("input data \n" + input.print_str() + "\n");

  // create mesh
  log("creating mesh\n");
  lm::ReplicatedMesh mesh(*comm);
  long N = long(1. / mesh_size);
  if (mesh_file != "")
    mesh.read(input.d_mesh_file);
  else
    lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                          1., 0., 1., lm::HEX8);

  // create equation system
  log("creating equation system\n");
  lm::EquationSystems eq_sys(mesh);

  // create spatial field of hydraulic conductivity
  auto &p_field = eq_sys.add_system<lm::ExplicitSystem>("K");
  p_field.add_variable("k", lm::CONSTANT, lm::MONOMIAL);

  eq_sys.init();

  // setting up perfusion points
  log("creating random perfusion points\n");
  auto bbox = lm::MeshTools::create_bounding_box(mesh);
  auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
  auto l = (bbox.min() - bbox.max()).norm();
  eq_sys.parameters.set<lm::Point>("center") = xc;
  eq_sys.parameters.set<double>("length") = l;

  // create outlet point data
  std::vector<lm::Point> pts;
  std::vector<double> radii;
  std::vector<double> flow;
  std::vector<double> vol;
  set_perfusion_pts(out_dir, pts, radii, hyd_cond, eq_sys, model);

  // create territory
  std::vector<long> partition;
  mc::estimate_flow_rate(radii, 3., flow);
  mc::estimate_vol_fraction(radii, flow, vol);
  mc::create_perfusion_territory(pts, vol,
                                 mesh, partition,
                                 200, hyd_cond, eq_sys,
                                 log, out_dir, 10);

  return 0;
}