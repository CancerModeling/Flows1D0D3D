////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "nifti_reader.hpp"

macrocirculation::NiftiReader::NiftiReader(std::string filename) {
  d_reader_p = vtkSmartPointer<vtkNIFTIImageReader>::New();
  d_reader_p->SetFileName(filename.c_str());
  d_reader_p->Update();

  d_img_p = vtkSmartPointer<vtkImageData>:: New();
  d_img_p->DeepCopy(d_reader_p->GetOutput());
}