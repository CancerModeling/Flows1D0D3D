////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_csv_writer.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

GraphCSVWriter::GraphCSVWriter(MPI_Comm comm,
                               std::string foldername,
                               std::string datasetname,
                               std::shared_ptr<GraphStorage> graph,
                               std::shared_ptr<DofMap> dof_map,
                               std::vector<std::string> component_names)
    : d_comm(comm),
      d_foldername(std::move(foldername)),
      d_datasetname(std::move(datasetname)),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_component_names(std::move(component_names)) {
  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto &edge = d_graph->get_edge(eid);
    const auto &local_dof_map = d_dof_map->get_local_dof_map(*edge);
    const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
    const auto h = length / local_dof_map.num_micro_edges();

    for (std::size_t component = 0; component < d_component_names.size(); component += 1) {
      std::fstream filecsv;
      filecsv.open(get_name(component, eid), std::ios::out);
      for (std::size_t micro_edge = 0; micro_edge < local_dof_map.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";
        filecsv << h * micro_edge << "," << h * (micro_edge + 1);
      }
      filecsv << std::endl;
    }
  }

  {
    std::fstream filecsv;
    filecsv.open(get_name_times(), std::ios::out);
    filecsv << "";
  }
}

std::string GraphCSVWriter::get_name(std::size_t component, std::size_t vessel) const {
  std::stringstream name;
  name << d_foldername << "/" << d_datasetname << "_" << d_component_names.at(component) << "_vessel" << std::setfill('0') << std::setw(5) << vessel
       << ".csv";
  return name.str();
}

std::string GraphCSVWriter::get_name_times() const {
  return d_foldername + "/" + d_datasetname + "_times.csv";
}

void GraphCSVWriter::write(double time, const std::vector<double> &data) const {
  // write time step information
  {
    std::fstream filecsv;
    filecsv.open(get_name_times(), std::ios::app);
    filecsv << time << std::endl;
  }

  // write vessel data
  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto &edge = d_graph->get_edge(eid);
    const auto &local_dof_map = d_dof_map->get_local_dof_map(*edge);
    const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
    const auto h = length / local_dof_map.num_micro_edges();

    FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions(), 0);
    std::vector<double> local_data(local_dof_map.num_basis_functions(), 0);

    for (std::size_t component = 0; component < d_component_names.size(); component += 1) {
      std::fstream filecsv;
      filecsv.open(get_name(component, eid), std::ios::app);
      for (std::size_t micro_edge = 0; micro_edge < local_dof_map.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";

        local_dof_map.dof_indices(micro_edge, component, dof_indices);
        extract_dof(dof_indices, data, local_data);
        auto boundary_values = fe.evaluate_dof_at_boundary_points(local_data);

        filecsv << boundary_values.left << "," << boundary_values.right;
      }
      filecsv << std::endl;
    }
  }
}

} // namespace macrocirculation