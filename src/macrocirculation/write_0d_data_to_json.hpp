////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_WRITE_0D_DATA_TO_JSON_HPP
#define TUMORMODELS_WRITE_0D_DATA_TO_JSON_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace macrocirculation {

// forward declarations:
class Values0DModel;
class GraphStorage;

void write_0d_data_to_json(const std::string &output_path,
                           const std::shared_ptr<GraphStorage> &graph,
                           const std::vector<double> &t,
                           const std::map<size_t, std::vector<Values0DModel>> &input);

} // namespace macrocirculation

#endif //TUMORMODELS_WRITE_0D_DATA_TO_JSON_HPP
