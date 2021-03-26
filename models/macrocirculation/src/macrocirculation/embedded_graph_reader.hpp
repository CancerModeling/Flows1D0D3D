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
    : d_elastic_modulus(0.8)
    , d_wall_width(0.034)
    , d_rho(1.028)
    , d_poisson_ratio(0.5)
  {}

  /*! @brief Appends the data from the input file to the given graph storage. */
  void append(const std::string &filepath, GraphStorage & graph) const;

  /*! @brief Sets some guesses for parameter not in the file. */
  void set_parameter(double poisson_ratio, double wall_width, double nu, double rho);

private:
  double d_elastic_modulus;
  double d_wall_width;
  double d_poisson_ratio;
  double d_rho;
};

} // namespace macrocirculation

#endif //TUMORMODELS_EMBEDDED_GRAPH_READER_HPP
