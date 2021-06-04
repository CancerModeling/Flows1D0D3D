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

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredPoints.h>
#include <vtkImageData.h>
#include <vtkDataObject.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedIntArray.h>

namespace mc = macrocirculation;

/*! @brief Given 3d vector element id, compute 1d element id for image data */
int index_3d_1d(int I[3], std::vector<int> dim) {

  return I[0] + I[1] * dim[0] + I[2] * dim[0] * dim[1];
}


int main(int argc, char *argv[]) {

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value(""))                                       //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/test.vtk"))                            //
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

  auto filename = args["mesh-file"].as<std::string>();

  // read file
  vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  reader->Print(std::cout);

  // get structured point data
  vtkSmartPointer<vtkStructuredPoints> data = vtkSmartPointer<vtkStructuredPoints>::New();
  data->DeepCopy(reader->GetOutput());
  data->Print(std::cout);

  int dim[3];
  data->GetDimensions(dim);
  std::cout << "data dimension = (" << dim[0] << ", " << dim[1] << ", " << dim[2] << ")\n";

  // get point data
  auto pt_data = data->GetPointData();
  pt_data->Print(std::cout);
  std::vector<std::string> fields;
  for (size_t i=0; i<pt_data->GetNumberOfArrays(); i++)
    fields.push_back(pt_data->GetArrayName(i));
  auto arr = pt_data->GetArray(fields[0].c_str());

  // read values
  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory

  {
    int ii[3] = {2, 14, 7};
    auto i1 = mc::index_3d_1d(ii, dim);
    auto i2 = data->ComputePointId(ii);
    arr->GetTuples(i2, i2, data_a);
    double a = data_a->GetValue(0);
    std::cout << "i (code) = " << i1 << ", i (vtk) = " << i2 << ", val 1 = " << a << "\n";

    ii[0] = 10;
    i1 = mc::index_3d_1d(ii, dim);
    i2 = data->ComputePointId(ii);
    arr->GetTuples(i2, i2, data_a);
    a = data_a->GetValue(0);
    std::cout << "i (code) = " << i1 << ", i (vtk) = " << i2 << ", val 2 = " << a << "\n";

    ii[0] = 14;
    ii[1] = 1;
    i1 = mc::index_3d_1d(ii, dim);
    i2 = data->ComputePointId(ii);
    arr->GetTuples(i2, i2, data_a);
    a = data_a->GetValue(0);
    std::cout << "i (code) = " << i1 << ", i (vtk) = " << i2 << ", val 3 = " << a << "\n";
  }

  return 0;
}