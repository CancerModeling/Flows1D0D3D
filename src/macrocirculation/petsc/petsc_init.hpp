////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "petsc.h"

#ifndef TUMORMODELS_PETSC_INIT_HPP
#define TUMORMODELS_PETSC_INIT_HPP

namespace macrocirculation {

class PetscInit {
public:
  PetscInit(int argc, char** argv){
    CHKERRABORT(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));
  }

  ~PetscInit(){
    CHKERRABORT(PETSC_COMM_WORLD, PetscFinalize());
  }

  PetscInit(const PetscInit&) = delete;
  PetscInit& operator=(const PetscInit&) = delete;
};

}


#endif //TUMORMODELS_PETSC_INIT_HPP
