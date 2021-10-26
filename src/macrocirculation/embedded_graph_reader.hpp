////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EMBEDDED_GRAPH_READER_HPP
#define TUMORMODELS_EMBEDDED_GRAPH_READER_HPP

#include <string>
#include <vector>

namespace macrocirculation {

class GraphStorage;

/*! @brief Reads the json meshes and appends them to graphs. */
class EmbeddedGraphReader {
public:
  EmbeddedGraphReader()
    : d_rho(1.028e-3)
  {}

  /*! @brief Appends the data from the input file to the given graph storage. */
  void append(const std::string &filepath, GraphStorage & graph) const;

  /*! @brief Overwrites the boundary data for the vertices in the given file. */
  void set_boundary_data(const std::string &filepath, GraphStorage & graph) const;

private:
  double d_rho;
};

struct InputPressuresResults {
  std::string name;
  std::vector< double > t;
  std::vector< double > p;
  bool periodic;
};

std::vector<InputPressuresResults> read_input_pressures(const std::string& filepath);

} // namespace macrocirculation

#endif //TUMORMODELS_EMBEDDED_GRAPH_READER_HPP
