////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "nifti_reader.hpp"

#include "vtkType.h"
#include "vtkPointData.h"
#include "vtkAbstractArray.h"
#include "vtkDoubleArray.h"
#include "vtkNIFTIImageHeader.h"

macrocirculation::NiftiReader::NiftiReader(std::string filename) {
  d_reader_p = vtkSmartPointer<vtkNIFTIImageReader>::New();
  d_reader_p->SetFileName(filename.c_str());
  d_reader_p->Update();

  d_img_p = vtkSmartPointer<vtkImageData>:: New();
  d_img_p->DeepCopy(d_reader_p->GetOutput());
}

std::string macrocirculation::NiftiReader::print_str() {
  std::ostringstream oss;
  oss << "Nifti header\n";
  d_reader_p->GetNIFTIHeader()->Print(oss);
  oss << "\n\nNifti data\n";
  d_img_p->Print(oss);
  oss << "\n\nPoint fields in the file\n";
  auto pt_data = d_img_p->GetPointData();
  oss << "number of fields = " << pt_data->GetNumberOfArrays() << "\n";
  for (size_t i=0; i<pt_data->GetNumberOfArrays(); i++)
    oss << i << " field name = " << pt_data->GetArrayName(i) << "\n";

  return oss.str();
}

std::vector<std::string> macrocirculation::NiftiReader::get_point_fields() {
  auto pt_data = d_img_p->GetPointData();
  std::vector<std::string> fields;
  for (size_t i=0; i<pt_data->GetNumberOfArrays(); i++)
    fields.push_back(pt_data->GetArrayName(i));
  return fields;
}


std::vector<int> macrocirculation::NiftiReader::get_data_dimension() {
  int dim[3];
  d_img_p->GetDimensions(dim);
  return {dim[0], dim[1], dim[2]};
}

void macrocirculation::NiftiReader::read_point_data(std::string field_name, std::vector<double> *data) {
  auto pt_data = d_img_p->GetPointData();
  if (!pt_data->HasArray(field_name.c_str())) {
    std::cerr << "Image file does not have point field with name = " << field_name << ". Exiting.\n";
    exit(EXIT_FAILURE);
  }
  auto arr = pt_data->GetArray(field_name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory

  // read
  data->resize(arr->GetNumberOfTuples());
  for (size_t i=0; i<arr->GetNumberOfTuples(); i++) {
    arr->GetTuples(i, i, data_a);
    (*data)[i] = data_a->GetValue(0);
  }
}

void macrocirculation::NiftiReader::read_point_data(std::string field_name, std::vector<std::vector<std::vector<double>>> *data) {
  auto pt_data = d_img_p->GetPointData();
  if (!pt_data->HasArray(field_name.c_str())) {
    std::cerr << "Image file does not have point field with name = " << field_name << ". Exiting.\n";
    exit(EXIT_FAILURE);
  }
  auto arr = pt_data->GetArray(field_name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory

  auto dim = get_data_dimension();
  (*data).resize(dim[0]);
  for (int i=0; i<dim[0]; i++) {
    (*data)[i].resize(dim[1]);
    for (int j=0; j<dim[1]; j++) {
      (*data)[i][j].resize(dim[2]);
      for (int k=0; k<dim[2]; k++) {
        int ijk = k + j * dim[2] + i * dim[2] * dim[1];
        arr->GetTuples(ijk, ijk, data_a);

        (*data)[i][j][k] = data_a->GetValue(0);
      }
    }
  }
}