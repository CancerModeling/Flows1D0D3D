////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP
#define TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP

#include <mpi.h>
#include <vector>
#include <functional>

namespace macrocirculation {

class GraphStorage;
class DofMap;
class Point;
class Edge;

void linear_interpolate_points(const Point &left,
                               const Point &right,
                               std::size_t num_micro_edges,
                               std::vector<Point> &points);

void interpolate_to_vertices(MPI_Comm comm,
                             const GraphStorage &graph,
                             const DofMap &map,
                             std::size_t component,
                             const std::vector<double> &dof_vector,
                             std::vector<Point> &points,
                             std::vector<double> &interpolated);

void add_discontinuous_points( const std::vector< Point >& embedded_points, std::vector< Point >& points );

void fill_with_vessel_id(const MPI_Comm comm, const GraphStorage &graph, std::vector<Point> &points, std::vector<double> &interpolated);

void fill_with_vessel_A0(const MPI_Comm comm, const GraphStorage &graph, std::vector<Point> &points, std::vector<double> &interpolated);

void fill_with_radius(const MPI_Comm comm, const GraphStorage &graph, std::vector<Point> &points, std::vector<double> &interpolated);

void fill_with_edge_parameter(const MPI_Comm comm, const GraphStorage &graph, std::function< double(const Edge&)> extractor, std::vector<Point> &points, std::vector<double> &interpolated);

} // namespace macrocirculation

#endif //TUMORMODELS_INTERPOLATE_TO_VERTICES_HPP
