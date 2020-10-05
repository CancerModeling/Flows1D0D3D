////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "AvaLima/model.hpp"
#include "AngLima/model.hpp"
#include "AvaMech/model.hpp"
#include "NetTum/model.hpp"
#include <libmesh/getpot.h>
#include <iostream>

int main(int argc, char *argv[]) {

  // Open file with model setup
  std::string filename = "input.in";

  // Open file with model setup
  GetPot input_file(filename);

  // find which model to run
  std::string model_name = input_file("model_name", "AvaLima");

  if (model_name == "AvaLima") {

    // Run the model
    std::vector<double> out_rd;
    new avalima::Model(argc, argv, out_rd, filename);

    // Print solution to screen
    for (auto i : out_rd)
      std::cout << "Solution: Tumor area = " << i << std::endl;

  } else if (model_name == "AngLima") {

    // Run the model
    std::vector<double> out_rd;
    new anglima::Model(argc, argv, out_rd, filename);

    // Print solution to screen
    for (auto i : out_rd)
      std::cout << "Solution: Tumor area = " << i << std::endl;

  } else if (model_name == "AvaMech") {

    // Run the model
    std::vector<double> out_rd;
    new avamech::Model(argc, argv, out_rd, filename);

    // Print solution to screen
    for (auto i : out_rd)
      std::cout << "Solution: Tumor area = " << i << std::endl;
  } else if (model_name == "NetTum") {

    // Note: This one requires pointer to comm and therefore we have to init
    // libmesh and then call the constructor of model
    LibMeshInit init(argc, argv);

    // Run the model
    std::vector<double> out_rd;
    auto model = nettum::Model(argc, argv, out_rd, filename, &init.comm());

    // Print solution to screen
    for (auto i : out_rd)
      std::cout << "Solution: Tumor area = " << i << std::endl;
  }

  // End Code
  return 0;
}
