////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "libmesh/point.h"
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

namespace macrocirculation {

inline void add_points(const std::vector<libMesh::Point> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto points = vtkSmartPointer<vtkPoints>::New();
  for (int i=0; i<x.size(); i++)
    points->InsertNextPoint(x[i](0), x[i](1), x[i](2));
  data->SetPoints(points);
}

template <class T>
inline void add_array(std::string name, const std::vector<T> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1] = {0.};
  for (const auto &xi : x) {
    value[0] = double(xi);
    array->InsertNextTuple(value);
  }

  // write
  data->GetPointData()->AddArray(array);
}

inline void add_vec_array(std::string name, const std::vector<libMesh::Point> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(3);
  array->SetName(name.c_str());

  array->SetComponentName(0, "x");
  array->SetComponentName(1, "y");
  array->SetComponentName(2, "z");

  double value[3] = {0., 0., 0.};
  for (const auto &xi : x) {
    value[0] = xi(0);
    value[1] = xi(1);
    value[2] = xi(2);
    array->InsertNextTuple(value);
  }

  // write
  data->GetPointData()->AddArray(array);
}

class VTKWriter {
 public:
  VTKWriter(std::string name) {
    d_w_p = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    d_d_p = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (!name.empty())
      set_filename(name);
  }

  void set_filename(std::string name) {d_w_p->SetFileName(name.c_str()); }

  void write() {
    d_w_p->SetInputData(d_d_p);
    d_w_p->SetDataModeToAppended();
    d_w_p->EncodeAppendedDataOn();
    d_w_p->Write();
  }

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> d_w_p;
  vtkSmartPointer<vtkUnstructuredGrid> d_d_p;
};


} // namespace macrocirculation

#endif
