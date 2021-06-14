////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_PVD_WRITER_HPP
#define TUMORMODELS_GRAPH_PVD_WRITER_HPP

#include <mpi.h>
#include <ostream>
#include <string>
#include <vector>

namespace macrocirculation {

// forward declarations
class Point;

/*! @brief Graph agnostic data to write serialized graph data to pvd files in parallel.
 *
 * Design considerations:
 * Q: Why does the class not know the graph?
 * A: By just providing the points we can only write to disk some part of the graph.
 *    This is relevant for us, since we are mixing a linear and a nonlinear model on it, which have different QoIs.
 *    Thus only the vertex interpolation routines have to know which parts to serialize and not the pvd writer.
 */
class GraphPVDWriter {
public:
  /*! @brief Constructs a pvd writer, which writes into folder_name  pvd, vtp and pvtp files with the given dataset_name. */
  GraphPVDWriter(MPI_Comm comm, std::string folder_name, std::string dataset_name);

  /*! @brief Sets the points of the graph we want to write next.
   *         The point data is ordered, such that the ith vessel segment is given by (points[2*i], points[2*i+1]).
   */
  void set_points(std::vector<Point> points);

  /*! @brief Adds vertex data to save later.
   *         The data has to match the serialized point data given in set_points.
   *         Thus depending on the connectivity every vertex gets saved several times.
   *         This makes sense, since our data is discontinuous anyway.
   */
  void add_vertex_data(const std::string &name, std::vector<double> data);

  /*! @brief Writes the data for the current time. */
  void write(double time);

private:
  const MPI_Comm d_comm;
  const std::string d_folder_name;
  const std::string d_dataset_name;

  std::vector<double> d_times;

  std::vector<Point> d_points;

  struct NamedField {
    NamedField(std::string n, std::vector<double> v) : name(std::move(n)), values(std::move(v)) {}

    std::string name;
    std::vector<double> values;
  };

  // Data-format: vertex_data.values[ edge_index * 2 + local_vertex_idx ]
  std::vector<NamedField> d_vertex_data;

private:
  void write_vtp(std::ostream &out) const;

  void write_pvtp(std::ostream &out) const;

  void write_pvd(std::ostream &out) const;

  std::string get_vtp_filename(std::size_t time_step, int rank) const;

  std::string get_pvtp_filename(std::size_t time_step) const;

  std::string get_pvd_filename() const;
};

} // namespace macrocirculation

#endif //TUMORMODELS_GRAPH_PVD_WRITER_HPP
