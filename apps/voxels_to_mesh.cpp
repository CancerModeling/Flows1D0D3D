////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cxxopts.hpp>
#include <fmt/format.h>
#include <memory>

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"
#include "macrocirculation/nifti_reader.hpp"
#include "macrocirculation/vtk_io_libmesh.hpp"

#include "macrocirculation/vtk_reader.hpp"

#include <zlib.h>

namespace mc = macrocirculation;

namespace darcy3d {

// read vascular domain from nifti file and create hydraulic parameter with one value in
// vascular domain and other outside the vascular domain
void create_heterogeneous_conductivity(std::string out_dir, std::string nifti_filename, double voxel_size, lm::MeshBase &mesh, lm::ExplicitSystem &hyd_cond, lm::EquationSystems &eq_sys) {

  // read file
  auto nifti = mc::NiftiReader(nifti_filename);
  auto img_dim = nifti.get_data_dimension();
  auto img_fields = nifti.get_point_fields();
  std::vector<double> img_data;
  nifti.read_point_data(img_fields[0], &img_data);
  std::vector<std::vector<std::vector<double>>> img_data_grid;
  nifti.read_point_data(img_fields[0], &img_data_grid);

  std::vector<unsigned int> dof_indices;
  for (const auto &elem : mesh.active_local_element_ptr_range()) {
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);
    auto x = elem->centroid() / voxel_size; // + lm::Point(0.5, 0.5, 0.5);

    // find voxel element (1d vector representation of the image data)
    int i = mc::locate_voxel_1d({x(0), x(1), x(2)}, img_dim);
    auto a = img_data[i];

    auto ii = mc::locate_voxel_3d({x(0), x(1), x(2)}, img_dim);
    auto xv = lm::Point(ii[0], ii[1], ii[2]) * voxel_size;

    // set parameter
    if (a > 0.5)
      hyd_cond.solution->set(dof_indices[0], 1);
    else
      hyd_cond.solution->set(dof_indices[0], 0.1);
  }
  hyd_cond.solution->close();
  hyd_cond.update();
}

/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
std::string compress_string(const std::string& str,
                            int compressionlevel = Z_BEST_COMPRESSION)
{
  z_stream zs;                        // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, compressionlevel) != Z_OK)
    throw(std::runtime_error("deflateInit failed while compressing."));

  zs.next_in = (Bytef*)str.data();
  zs.avail_in = str.size();           // set the z_stream's input

  int ret;
  char outbuffer[32768];
  std::string outstring;

  // retrieve the compressed bytes blockwise
  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = deflate(&zs, Z_FINISH);

    if (outstring.size() < zs.total_out) {
      // append the block to the output string
      outstring.append(outbuffer,
                       zs.total_out - outstring.size());
    }
  } while (ret == Z_OK);

  deflateEnd(&zs);

  if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return outstring;
}

/** Decompress an STL string using zlib and return the original data. */
std::string decompress_string(const std::string& str)
{
  z_stream zs;                        // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (inflateInit(&zs) != Z_OK)
    throw(std::runtime_error("inflateInit failed while decompressing."));

  zs.next_in = (Bytef*)str.data();
  zs.avail_in = str.size();

  int ret;
  char outbuffer[32768];
  std::string outstring;

  // get the decompressed bytes blockwise using repeated calls to inflate
  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = inflate(&zs, 0);

    if (outstring.size() < zs.total_out) {
      outstring.append(outbuffer,
                       zs.total_out - outstring.size());
    }

  } while (ret == Z_OK);

  inflateEnd(&zs);

  if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") "
        << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return outstring;
}
} // namespace darcy3d

int main(int argc, char *argv[]) {
  auto sim_begin = std::chrono::steady_clock::now();

  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/meshes/darcy_test_tissue_mesh.e"))                 //
    ("nifti-file", "nifti file containing parameter values", cxxopts::value<std::string>()->default_value("data/meshes/darcy_test_vascular_domain.nii.gz")) //
    ("voxel-size", "voxel size used in creating the tissue mesh", cxxopts::value<double>()->default_value("1."))                                            //
    ("output-dir", "directory for the output", cxxopts::value<std::string>()->default_value("./output_voxels_to_mesh/"))                                    //
    ("output-file", "filename for the output", cxxopts::value<std::string>()->default_value("data"))                                                        //
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

  auto mesh_file = args["mesh-file"].as<std::string>();
  auto out_dir = args["output-dir"].as<std::string>();
  auto out_file = out_dir + args["output-file"].as<std::string>();
  auto nifti_file = args["nifti-file"].as<std::string>();
  auto voxel_size = args["voxel-size"].as<double>();

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  // create mesh
  log("creating mesh\n");
  lm::ReplicatedMesh mesh(*comm);


  if (true) {
    mesh.read(mesh_file);

    log("creating equation system\n");
    lm::EquationSystems eq_sys(mesh);
    auto &hyd_cond = eq_sys.add_system<lm::ExplicitSystem>("K");
    hyd_cond.add_variable("k", lm::CONSTANT, lm::MONOMIAL);
    eq_sys.init();

    // create heterogeneous property field
    log("setting up K field\n");
    auto bbox = lm::MeshTools::create_bounding_box(mesh);
    auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
    auto l = (bbox.min() - bbox.max()).norm();
    eq_sys.parameters.set<lm::Point>("center") = xc;
    eq_sys.parameters.set<double>("length") = l;

    darcy3d::create_heterogeneous_conductivity(out_dir, nifti_file, voxel_size, mesh, hyd_cond, eq_sys);

    // write
    log("writing to file\n");
    std::set< std::string > str_set = {"k"};
    mc::VTKIO(mesh).write_equation_systems(out_dir + "output_0.pvtu", eq_sys, &str_set);
    {
      auto temp_fname = out_dir + "output_0.e";
      lm::ExodusII_IO exo(mesh);
      exo.write_equation_systems(temp_fname.c_str(), eq_sys);
      exo.write_timestep(temp_fname.c_str(), eq_sys, 1, 0.);
    }
  }

  // test libmesh vtu reader
  auto fname = out_dir + "output_0_0.vtu";
  mesh.read(fname);
  mesh.print_info(lm::out);

  // test custom vtu reader
  auto reader = mc::VTKReader(fname);
  reader.read();
  auto p_fields = mc::read_vtu_point_fields(reader.d_d_p);
  auto c_fields = mc::read_vtu_point_fields(reader.d_d_p);

  log(fmt::format("\n\nPoint fields = {}\n", mc::print_str(p_fields)));
  log(fmt::format("\nCell fields = {}\n", mc::print_str(c_fields)));

  std::vector<double> p_k;
  mc::read_point_array("k", p_k, reader.d_d_p);
  log(fmt::format("Size of point k field = {}\n", p_k.size()));
  std::vector<double> c_k;
  mc::read_point_array("k", c_k, reader.d_d_p);
  log(fmt::format("Size of cell k field = {}\n", c_k.size()));

  if (p_k.size() == c_k.size()) {
    std::vector<double> pc_k;
    for (size_t i = 0; i < p_k.size(); i++)
      pc_k.push_back(p_k[i] - c_k[i]);

    log(fmt::format("max err = {}, min err = {}\n", mc::avg(pc_k), mc::std(pc_k)));

    // find if there is a element with value in between 0.1 and 1
    bool found_intermed_elems = false;
    for (auto a : p_k) {
      if (a > 0.1 + 1.e-5 and a < 1. - 1.e-5)
        log(fmt::format("a = {}\n", a));
    }
  }

  return 0;
}