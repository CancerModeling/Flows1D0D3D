////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_EMBEDDED_GRAPH_READER_HPP
#define TUMORMODELS_EMBEDDED_GRAPH_READER_HPP

#include <string>

namespace macrocirculation {

class GraphStorage;

class EmbeddedGraphReader {
public:
  EmbeddedGraphReader()
    : d_rho(1.028e-3)
  {}

  /*! @brief Appends the data from the input file to the given graph storage. */
  void append(const std::string &filepath, GraphStorage & graph) const;

  /*! @brief Sets some guesses for parameter not in the file. */
  void set_parameter(double poisson_ratio, double wall_width, double nu, double rho);

private:
  double d_rho;
};

} // namespace macrocirculation

#endif //TUMORMODELS_EMBEDDED_GRAPH_READER_HPP
