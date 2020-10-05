////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "model.hpp"
#include <libmesh/getpot.h>
#include <iostream>

int main(int argc, char *argv[]) {

  // Open file with model setup
  std::string filename = "input.in";

  // Open file with model setup
  GetPot input_file(filename);

  // find which model to run
  std::string model_name = input_file("model_name", "TwoSpecies");

  if (model_name == "TwoSpecies") {

    // Note: This one requires pointer to comm and therefore we have to init
    // libmesh and then call the constructor of model
    LibMeshInit init(argc, argv);

    // Run the model
    twosp::model_setup_run(argc, argv, filename, &init.comm());

  } else {
    std::cout << "Model name = " << model_name
              << " can not be run using TwoSpecies driver.\n ";
  }

  // End Code
  return 0;
}
