////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NIFTI_READER_HPP
#define NIFTI_READER_HPP

#include <vector>
#include "vtkSmartPointer.h"
#include "vtkNIFTIImageReader.h"
#include "vtkImageData.h"

namespace macrocirculation {

/*! @brief Simple API for reading nifti files using vtk libraries */
class NiftiReader{

public:
  /*! @brief Constructs a pvd writer, which writes into folder_name  pvd, vtp and pvtp files with the given dataset_name. */
  NiftiReader(std::string filename);

  /*! @brief get data information */
  std::string print_str();

  /*! @brief Get fields stored in the data. */
  std::vector<std::string> get_point_fields();

  /*! @brief Get fields stored in the data. */
  std::vector<int> get_data_dimension();

  /*! @brief Read data into vector. */
  void read_point_data(std::string field_name, std::vector<double> *data);

  /*! @brief Read data into 3 dimensional array. */
  void read_point_data(std::string field_name, std::vector<std::vector<std::vector<double>>> *data);

  /*! @brief vtk reader */
  vtkSmartPointer<vtkNIFTIImageReader> d_reader_p;

  /*! @brief image data */
  vtkSmartPointer<vtkImageData> d_img_p;
};


} // namespace macrocirculation

#endif //NIFTI_READER_HPP
