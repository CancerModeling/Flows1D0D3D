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

#include "macrocirculation/nifti_reader.hpp"

namespace mc = macrocirculation;

int main(int argc, char *argv[]) {
  cxxopts::Options options(argv[0], "Test nifti reader");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value("data/test.nii.gz")) //
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
  std::cout << "Print nifti header file\n";
  read_img.d_reader_p->GetNIFTIHeader()->Print(std::cout);
  std::cout << "Print nifti data\n";
  read_img.d_reader_p->GetOutput()->Print(std::cout);

  auto pt_data = read_img.d_img_p->GetPointData();
  auto arr = pt_data->GetArray("NIFTI");
  std::cout << "Print NIFTI data in nifti file\n";
  arr->Print(std::cout);
  // read and print data
  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory
  for (size_t i=0; i<arr->GetNumberOfTuples(); i++) {
    arr->GetTuples(i, i, data_a);
    std::cout << i << ", " << data_a->GetValue(0) << "\n";
  }

  return 0;
}