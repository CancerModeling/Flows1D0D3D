////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_data_writer.hpp"

#include "graph_storage.hpp"
#include "dof_map_network.hpp"

#include <fstream>

namespace macrocirculation {

void GraphDataWriter::add_midpoint_data(const std::string &name, std::vector<double> data) {
  midpoint_data.push_back(std::move(NamedField(name, std::move(data))));
}

std::string GraphDataWriter::get_path(const std::string &filename, std::size_t time_idx) const {
  std::string path = filename + "_";
  path += std::to_string(time_idx);
  path += ".vtk";
  return path;
}

void GraphDataWriter::write_vtk(const std::string &filename, GraphStorage &storage, std::size_t time_idx) const {
  const auto path = get_path(filename, time_idx);

  const std::size_t number_of_points = storage.num_edges() * 2;
  const std::size_t number_of_segments = storage.num_edges();

  std::fstream filevtk;
  filevtk.open(path, std::ios::out);
  filevtk << "# vtk DataFile Version 2.0" << std::endl;
  filevtk << "Network Nutrient Transport" << std::endl;
  filevtk << "ASCII" << std::endl;
  filevtk << "DATASET POLYDATA" << std::endl;

  filevtk << "POINTS " << number_of_points << " float" << std::endl;
  // the part for our point data
  for (auto e_id : storage.get_edge_ids()) {
    auto edge = storage.get_edge(e_id);
    const auto &p0 = edge->get_coordinate_v0();
    const auto &p1 = edge->get_coordinate_v1();

    // write the two points into our file
    filevtk << p0(0) << " " << p0(1) << " " << p0(2) << std::endl;
    filevtk << p1(0) << " " << p1(1) << " " << p1(2) << std::endl;
  }
  filevtk << " " << std::endl;

  filevtk << "LINES " << number_of_segments << " " << 3 * number_of_segments << std::endl;
  for (std::size_t idx = 0; idx < number_of_segments; idx += 1)
    filevtk << "2 " << 2 * idx << " " << 2 * idx + 1 << std::endl;
  filevtk << " " << std::endl;

  filevtk << "CELL_DATA " << number_of_segments << std::endl;
  for (std::size_t idx = 0; idx < midpoint_data.size(); idx += 1) {
    filevtk << "SCALARS " << midpoint_data[idx].name << " float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    for (auto e_id : storage.get_edge_ids()) {
      auto edge = storage.get_edge(e_id);
      filevtk << midpoint_data[idx].values[edge->get_id()] << std::endl;
    }
  }
  filevtk << " " << std::endl;

  filevtk << "POINT_DATA " << number_of_points << std::endl;
  for (std::size_t idx = 0; idx < vertex_data.size(); idx += 1) {
    filevtk << "SCALARS " << vertex_data[idx].name << " float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    for (auto e_id : storage.get_edge_ids()) {
      filevtk << vertex_data[idx].values[2*e_id + 0] << std::endl;
      filevtk << vertex_data[idx].values[2*e_id + 1] << std::endl;
    }
  }
}

void GraphDataWriter::add_vertex_data(const std::string &name, std::vector<double> data)
{
  vertex_data.push_back(std::move(NamedField(name, std::move(data))));
}

} // namespace macrocirculation