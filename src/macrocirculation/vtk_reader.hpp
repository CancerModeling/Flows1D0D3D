////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef VTK_READER_HPP
#define VTK_READER_HPP

#include "libmesh/point.h"
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

namespace macrocirculation {

    /*!
   * @brief Checks if file has needed data
   * @param name Name of the variable
   * @return True If file has variable otherwise false
   */
    bool vtu_has_point_data(const std::string &name) {};

    /*!
   * @brief Checks if file has needed data
   * @param name Name of the variable
   * @return True If file has variable otherwise false
   */
    bool vtu_has_cell_data(const std::string &data_tag){};

    /*!
     * @brief Creates list of available data
     * @return List List of name of point data
     */
    std::vector<std::string> read_vtu_point_fields(){};

    /*!
     * @brief Creates list of available data
     * @return List List of name of cell data
     */
    std::vector<std::string> read_vtu_cell_fields(){};

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
