////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <utility>

#include "communication/mpi.hpp"
#include "graph_pvd_writer.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

GraphPVDWriter::GraphPVDWriter(MPI_Comm comm, std::string folder_name, std::string dataset_name)
    : d_comm(comm),
      d_folder_name(std::move(folder_name)),
      d_dataset_name(std::move(dataset_name)),
      d_times(),
      d_points(),
      d_vertex_data() {}

void GraphPVDWriter::set_points(std::vector<Point> points) {
  if (!d_vertex_data.empty())
    throw std::runtime_error("points data must be set before vertex data");
  d_points = std::move(points);
}

void GraphPVDWriter::add_vertex_data(const std::string &name, std::vector<double> data) {
  if (d_points.size() != data.size())
    throw std::runtime_error("vector size " + std::to_string(data.size()) +
                             " of " + name + " does not match number of points " + std::to_string(d_points.size()));

  d_vertex_data.push_back(std::move(NamedField(name, std::move(data))));
}

void GraphPVDWriter::write(double time) {
  if (!d_times.empty() && time < d_times.back())
    throw std::runtime_error("inserting values of the past in graph writer.");

  d_times.push_back(time);

  // write vtp file for current rank
  {
    std::fstream vtp_file(d_folder_name + "/" + get_vtp_filename(d_times.size() - 1, mpi::rank(d_comm)), std::ios::out);
    write_vtp(vtp_file);
  }

  // write other files if root
  if (mpi::rank(d_comm) == 0) {
    // write pvtp file
    {
      std::fstream pvtp_file(d_folder_name + "/" + get_pvtp_filename(d_times.size() - 1), std::ios::out);
      write_pvtp(pvtp_file);
    }

    // write pvd file
    {
      std::fstream pvd_file(d_folder_name + "/" + get_pvd_filename(), std::ios::out);
      write_pvd(pvd_file);
    }
  }

  // clear all vertex data, though we cache the point data for later.
  d_vertex_data.clear();
}

void GraphPVDWriter::write_vtp(std::ostream &out) const {
  // check that the data size and point size match
  for (auto &nf : d_vertex_data)
    if (nf.values.size() != d_points.size())
      throw std::runtime_error("vector size " + std::to_string(nf.values.size()) +
                               " of " + nf.name + " does not match number of points " + std::to_string(d_points.size()));

  const std::size_t num_points = d_points.size();

  if (num_points % 2 != 0)
    throw std::runtime_error("odd number of points cannot be divided into segments");

  const std::size_t num_lines = num_points / 2;

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  out << "<PolyData>\n";
  out << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfLines=\"" << num_lines << "\">\n";

  // write the point coordinates
  out << "<Points>\n";
  out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : d_points)
    out << p.x << " " << p.y << " " << p.z << " ";
  out << "\n";
  out << "</DataArray>\n";
  out << "</Points>\n";

  // connect the points with lines
  out << "<Lines>\n";
  // the connectivity
  out << "<DataArray type=\"Int32\" Name=\"connectivity\">\n";
  for (std::size_t idx = 0; idx < num_points; idx += 1)
    out << std::to_string(idx) << " ";
  out << "\n";
  out << "</DataArray>\n";
  // the end positions of all the cells (here line segments)
  out << "<DataArray type=\"Int32\" Name=\"offsets\">\n";
  for (std::size_t idx = 0; idx < num_lines; idx += 1)
    out << std::to_string(2 * (idx + 1)) << " ";
  out << "\n";
  out << "</DataArray>\n";
  out << "</Lines>\n";

  // write the scalar data
  if (!d_vertex_data.empty()) {
    out << "<PointData Scalars=\"" << d_vertex_data[0].name << "\">\n";
    for (const auto &nf : d_vertex_data) {
      out << "<DataArray type=\"Float32\" format=\"ascii\" Name=\"" << nf.name << "\">\n";
      for (const auto &value : nf.values)
        out << value << " ";
      out << "\n";
      out << "</DataArray>\n";
    }
    out << "</PointData>\n";
  }

  out << "</Piece>\n";
  out << "</PolyData>\n";

  out << "</VTKFile>\n";
}

void GraphPVDWriter::write_pvtp(std::ostream &out) const {
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "<PPolyData GhostLevel=\"0\">\n";

  out << "<PPoints>\n";
  out << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\" />\n";
  out << "</PPoints>\n";

  if (!d_vertex_data.empty()) {
    out << "<PPointData Scalars=\"" << d_vertex_data[0].name << "\">\n";
    for (const auto &nf : d_vertex_data)
      out << "<PDataArray type=\"Float32\" Name=\"" << nf.name << "\" />\n";
    out << "</PPointData>\n";
  }

  for (std::size_t rank = 0; rank < static_cast< size_t > ( mpi::size(d_comm) ); rank += 1) {
    const auto filename = get_vtp_filename(d_times.size() - 1, rank);
    out << "<Piece Source=\"" << filename << "\" />\n";
  }

  out << "</PPolyData>\n";
  out << "</VTKFile>\n";
}

void GraphPVDWriter::write_pvd(std::ostream &out) const {
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "<Collection>\n";

  for (std::size_t time_index = 0; time_index < d_times.size(); time_index += 1) {
    const auto filename = get_pvtp_filename(time_index);
    out << "<DataSet timestep=\"" << d_times[time_index] << "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n";
  }

  out << "</Collection>\n";
  out << "</VTKFile>\n";
}

std::string GraphPVDWriter::get_vtp_filename(std::size_t time_step, int rank) const {
  return d_dataset_name + "_" + std::to_string(rank) + "_" + std::to_string(time_step) + ".vtp";
}

std::string GraphPVDWriter::get_pvtp_filename(std::size_t time_step) const {
  return d_dataset_name + "_" + std::to_string(time_step) + ".pvtp";
}

std::string GraphPVDWriter::get_pvd_filename() const {
  return d_dataset_name + ".pvd";
}

} // namespace macrocirculation
