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

#include "macrocirculation/utils.hpp"
#include "macrocirculation/nifti_reader.hpp"

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  cxxopts::Options options(argv[0], "Test nifti reader");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value("data/meshes/test.nii.gz")) //
    ("h,help", "print usage");
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if (!args.unmatched().empty()) {
    std::cout << "The following arguments were unmatched: " << std::endl;
    for (auto &it : args.unmatched())
      std::cout << " " << it;
  }

  auto filename = args["input-file"].as<std::string>();

  // read
  std::cout << "Reading nifti file = " << filename << "\n";
  auto read_img = mc::NiftiReader(filename);
  std::cout << "Print nifti file info\n";
  std::cout << read_img.print_str() << "\n";

  // read and print data
  std::vector<double> pt_data;
  auto field_names = read_img.get_point_fields();
  read_img.read_point_data(field_names[0], &pt_data);
  for (size_t i=0; i<pt_data.size(); i++) {
    std::cout << i << ", " << pt_data[i] << "\n";
  }

  // get image dimension
  auto dim = read_img.get_data_dimension();
  std::cout << "Dimension of the data = (" << dim[0] << ", " << dim[1] << ", " << dim[2] << ")\n";

  // check if grid wise data matches with 1d vector data
  std::vector<std::vector<std::vector<double>>> pt_data_grid;
  read_img.read_point_data(field_names[0], &pt_data_grid);
  size_t count = 0;
  for (int I=0; I<pt_data.size(); I++) {
    auto x = pt_data[I];
    auto I_3d = mc::index_1d_3d(I, dim);
    auto y = pt_data_grid[I_3d[0]][I_3d[1]][I_3d[2]];
    std::cout << "1d index = " << I
              << ", 3d index = (" << I_3d[0] << ", " << I_3d[1] << ", " << I_3d[2]
              << "), 1d value = " << x << ", 3d value = " << y << "\n";
    if (x != y)
      count++;
  }
  std::cout << "number of data not matching = " << count << "\n";

  return 0;
}