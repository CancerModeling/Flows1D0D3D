////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "boundary_condition_readers.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

#include "nonlinear_linear_coupling.hpp"

namespace macrocirculation {

// Expected json file structure
// {
//   "file_start": "<filename of first coupling file>",
//   "file_end": "<filename of second coupling file>",
//   "couplings": [{
//     "start": "<vertex-name>",
//     "end": "<vertex-name>"
//   }, ...]
// }
void read_coupling_conditions(NonlinearLinearCoupling & coupling, const std::string& filepath)
{
  using json = nlohmann::json;

  std::fstream file(filepath, std::ios::in);
  json j;
  file >> j;

  std::cout << "coupling between " << j["file_start"] << " and " << j["file_end"] << std::endl;

  auto couplings = j["couplings"];
  for (auto coupling_pair: couplings)
  {
    std::cout << "couples node " <<coupling_pair["start"] << " with node " << coupling_pair["end"] << std::endl;
    coupling.add_coupled_vertices( coupling_pair["start"], coupling_pair["end"] );
  }
}

}

