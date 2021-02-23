////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_DATA_WRITER_HPP
#define TUMORMODELS_GRAPH_DATA_WRITER_HPP

#include <string>
#include <vector>

namespace macrocirculation {

// forward declaration
class GraphStorage;

class GraphDataWriter {
public:
  void add_midpoint_data(const std::string &name, std::vector<double> data);

  /*! @brief Adds vertex data to save later.
   *         Note that the data format has to have the form
   *            vertex_data.values[ edge_index * 2 + local_vertex_idx ]
   *         and thus every vertex is saved several times, since our data is discontinuous.
   */
  void add_vertex_data(const std::string &name, std::vector<double> data);

  void write_vtk(const std::string &filename, GraphStorage &storage, std::size_t time_idx) const;

private:
  struct NamedField {
    NamedField(std::string n, std::vector<double> v) : name(std::move(n)), values(std::move(v)) {}

    std::string name;
    std::vector<double> values;
  };

  std::vector<NamedField> midpoint_data;

  // Data-format: vertex_data.values[ edge_index * 2 + local_vertex_idx ]
  std::vector<NamedField> vertex_data;

  std::string get_path(const std::string &filename, std::size_t time_idx) const;
};

} // namespace macrocirculation

#endif //TUMORMODELS_GRAPH_DATA_WRITER_HPP
