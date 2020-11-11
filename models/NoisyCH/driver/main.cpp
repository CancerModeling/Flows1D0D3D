////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/getpot.h"
#include "model.hpp"
#include <iostream>

int main(int argc, char *argv[]) {

  // Open file with model setup
  std::string filename = "input.in";

  // Note: This one requires pointer to comm and therefore we have to init
  // libmesh and then call the constructor of model
  LibMeshInit init(argc, argv);

  // Run the model
  noisych::model_setup_run(filename, &init.comm());

  // End Code
  return 0;
}
