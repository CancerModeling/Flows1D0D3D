////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef VTK_READER_HPP
#define VTK_READER_HPP

#include "libmesh/point.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

namespace macrocirculation {

/*!
 * @brief Checks if file has needed data
 * @param name Name of the variable
 * @param data VTU data
 * @return True If file has variable otherwise false
 */
bool vtu_has_point_data(const std::string &name, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto p_data = data->GetPointData();
  for (size_t i = 0; i < p_data->GetNumberOfArrays(); i++) {
    if (p_data->GetArrayName(i) == name)
      return true;
  }
  return false;
}

/*!
 * @brief Checks if file has needed data
 * @param name Name of the variable
 * @param data VTU data
 * @return True If file has variable otherwise false
 */
bool vtu_has_cell_data(const std::string &name, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto p_data = data->GetCellData();
  for (size_t i = 0; i < p_data->GetNumberOfArrays(); i++) {
    if (p_data->GetArrayName(i) == name)
      return true;
  }
  return false;
}

/*!
 * @brief Creates list of available data
 * @param data VTU data
 * @return List List of name of point data
 */
std::vector<std::string> read_vtu_point_fields(vtkSmartPointer<vtkUnstructuredGrid> &data) {
  std::vector<std::string> fields;
  auto p_data = data->GetPointData();
  for (size_t i = 0; i < p_data->GetNumberOfArrays(); i++)
    fields.push_back(p_data->GetArrayName(i));
  return fields;
}

/*!
 * @brief Creates list of available data
 * @param data VTU data
 * @return List List of name of cell data
 */
std::vector<std::string> read_vtu_cell_fields(vtkSmartPointer<vtkUnstructuredGrid> &data) {
  std::vector<std::string> fields;
  auto p_data = data->GetCellData();
  for (size_t i = 0; i < p_data->GetNumberOfArrays(); i++)
    fields.push_back(p_data->GetArrayName(i));
  return fields;
}

inline void read_points(std::vector<libMesh::Point> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {
  auto n = data->GetNumberOfPoints();
  x.clear();
  x.resize(n);
  double xi[3] = {0., 0., 0.};
  for (size_t i = 0; i < n; i++) {
    data->GetPoint(vtkIdType(i), xi);
    x[i] = libMesh::Point(xi[0], xi[1], xi[2]);
  }
}

template <class T>
inline bool read_point_array(std::string name, std::vector<T> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {

  auto p_data = data->GetPointData();
  if (p_data->HasArray(name.c_str()) == 0)
    return false;

  vtkDataArray *array = p_data->GetArray(name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory

  x.resize(array->GetNumberOfTuples());
  for (size_t i = 0; i < array->GetNumberOfTuples(); i++) {
    array->GetTuples(i, i, data_a);
    x[i] = data_a->GetValue(0);
  }

  return true;
}

template <class T>
inline bool read_cell_array(std::string name, std::vector<T> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {

  auto p_data = data->GetCellData();
  if (p_data->HasArray(name.c_str()) == 0)
    return false;

  vtkDataArray *array = p_data->GetArray(name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(1, 1);  // allocate memory

  x.resize(array->GetNumberOfTuples());
  for (size_t i = 0; i < array->GetNumberOfTuples(); i++) {
    array->GetTuples(i, i, data_a);
    x[i] = data_a->GetValue(0);
  }

  return true;
}

inline bool read_point_vec_array(std::string name, std::vector<libMesh::Point> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {

  auto p_data = data->GetPointData();
  if (p_data->HasArray(name.c_str()) == 0)
    return false;

  vtkDataArray *array = p_data->GetArray(name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(3, 1);  // allocate memory

  x.resize(array->GetNumberOfTuples());
  for (size_t i = 0; i < array->GetNumberOfTuples(); i++) {
    array->GetTuples(i, i, data_a);
    x[i] = libMesh::Point(data_a->GetValue(0), data_a->GetValue(1), data_a->GetValue(2));
  }

  return true;
}

inline bool read_cell_vec_array(std::string name, std::vector<libMesh::Point> &x, vtkSmartPointer<vtkUnstructuredGrid> &data) {

  auto p_data = data->GetCellData();
  if (p_data->HasArray(name.c_str()) == 0)
    return false;

  vtkDataArray *array = p_data->GetArray(name.c_str());

  auto data_a = vtkSmartPointer<vtkDoubleArray>::New();
  data_a->SetNumberOfComponents(1);
  data_a->Allocate(3, 1);  // allocate memory

  x.resize(array->GetNumberOfTuples());
  for (size_t i = 0; i < array->GetNumberOfTuples(); i++) {
    array->GetTuples(i, i, data_a);
    x[i] = libMesh::Point(data_a->GetValue(0), data_a->GetValue(1), data_a->GetValue(2));
  }

  return true;
}

class VTKReader {
public:
  VTKReader(std::string name) {
    d_r_p = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    d_d_p = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (!name.empty())
      set_filename(name);
  }

  void set_filename(std::string name) {
    d_r_p->SetFileName(name.c_str());
    d_r_p->Update();
  }

  void read() {
    d_d_p = d_r_p->GetOutput();
  }

  vtkSmartPointer<vtkXMLUnstructuredGridReader> d_r_p;
  vtkSmartPointer<vtkUnstructuredGrid> d_d_p;
};


} // namespace macrocirculation

#endif
