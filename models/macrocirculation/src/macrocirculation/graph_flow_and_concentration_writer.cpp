////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_flow_and_concentration_writer.hpp"

#include <cmath>
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

GraphFlowAndConcentrationWriter::GraphFlowAndConcentrationWriter(MPI_Comm comm,
                                                                 std::string foldername,
                                                                 std::string datasetname,
                                                                 std::shared_ptr<GraphStorage> graph,
                                                                 std::shared_ptr<DofMap> dof_map_flow,
                                                                 std::shared_ptr<DofMap> dof_map_transport)
    : d_comm(comm),
      d_foldername(std::move(foldername)),
      d_datasetname(std::move(datasetname)),
      d_graph(std::move(graph)),
      d_dof_map_flow(std::move(dof_map_flow)),
      d_dof_map_transport(std::move(dof_map_transport)) {
  const std::array<std::string, 2> flow_component_names = {"Q", "A"};

  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto &edge = d_graph->get_edge(eid);

    // prepare flow
    {
      const auto &local_dof_map_flow = d_dof_map_flow->get_local_dof_map(*edge);
      const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
      const auto h = length / local_dof_map_flow.num_micro_edges();

      for (std::size_t component = 0; component < 2; component += 1) {
        std::fstream filecsv;
        filecsv.open(get_name(flow_component_names[component], eid), std::ios::out);
        for (std::size_t micro_edge = 0; micro_edge < local_dof_map_flow.num_micro_edges(); micro_edge += 1) {
          if (micro_edge > 0)
            filecsv << ",";
          filecsv << h * micro_edge << "," << h * (micro_edge + 1);
        }
        filecsv << std::endl;
      }

      {
        std::fstream filecsv;
        filecsv.open(get_name("p", eid), std::ios::out);
        for (std::size_t micro_edge = 0; micro_edge < local_dof_map_flow.num_micro_edges(); micro_edge += 1) {
          if (micro_edge > 0)
            filecsv << ",";
          filecsv << h * micro_edge << "," << h * (micro_edge + 1);
        }
        filecsv << std::endl;
      }
    }

    // prepare concentration
    {
      const auto &local_dof_map_transport = d_dof_map_transport->get_local_dof_map(*edge);
      const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
      const auto h = length / local_dof_map_transport.num_micro_edges();

      std::fstream filecsv;
      filecsv.open(get_name("c", eid), std::ios::out);
      for (std::size_t micro_edge = 0; micro_edge < local_dof_map_transport.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";
        filecsv << h * micro_edge << "," << h * (micro_edge + 1);
      }
      filecsv << std::endl;
    }
  }

  // reset time
  if (mpi::rank(d_comm) == 0) {
    std::fstream filecsv;
    filecsv.open(get_name_times(), std::ios::out);
    filecsv << "";
  }
}

std::string GraphFlowAndConcentrationWriter::get_name(const std::string& name_quantity, std::size_t vessel) const {
  std::stringstream name;
  name << d_foldername << "/" << d_datasetname << "_" << name_quantity << "_vessel" << std::setfill('0') << std::setw(5) << vessel
       << ".csv";
  return name.str();
}

std::string GraphFlowAndConcentrationWriter::get_name_times() const {
  return d_foldername + "/" + d_datasetname + "_times.csv";
}

void GraphFlowAndConcentrationWriter::write(double time, const std::vector<double> &flow, const std::vector<double> &transport) const {
  // write time step information
  if (mpi::rank(d_comm) == 0) {
    std::fstream filecsv;
    filecsv.open(get_name_times(), std::ios::app);
    filecsv << time << std::endl;
  }
  const std::array<std::string, 2> flow_component_names = {"Q", "A"};

  // write vessel data
  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    {

      const auto &edge = d_graph->get_edge(eid);
      const auto &local_dof_map_flow = d_dof_map_flow->get_local_dof_map(*edge);
      const auto &local_dof_map_transport = d_dof_map_transport->get_local_dof_map(*edge);
      const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
      const auto h = length / local_dof_map_flow.num_micro_edges();

      FETypeNetwork fe(create_midpoint_rule(), local_dof_map_flow.num_basis_functions() - 1);
      fe.reinit(h);

      std::vector<std::size_t> dof_indices(local_dof_map_flow.num_basis_functions(), 0);
      std::vector<double> local_data(local_dof_map_flow.num_basis_functions(), 0);

      for (std::size_t component = 0; component < 2; component += 1) {
        std::fstream filecsv;
        filecsv.open(get_name(flow_component_names[component], eid), std::ios::app);
        for (std::size_t micro_edge = 0; micro_edge < local_dof_map_flow.num_micro_edges(); micro_edge += 1) {
          if (micro_edge > 0)
            filecsv << ",";

          local_dof_map_flow.dof_indices(micro_edge, component, dof_indices);
          extract_dof(dof_indices, flow, local_data);
          auto boundary_values = fe.evaluate_dof_at_boundary_points(local_data);

          filecsv << boundary_values.left << "," << boundary_values.right;
        }
        filecsv << std::endl;
      }

      // pressure
      {
        const auto A0 = edge->get_physical_data().A0;
        const auto G0 = edge->get_physical_data().G0;

        const auto pressure = [=](auto A) { return G0 * (std::sqrt(A / A0) - 1.0) / 1.33332; };

        std::fstream filecsv;
        filecsv.open(get_name("p", eid), std::ios::app);
        for (std::size_t micro_edge = 0; micro_edge < local_dof_map_flow.num_micro_edges(); micro_edge += 1) {
          if (micro_edge > 0)
            filecsv << ",";

          local_dof_map_flow.dof_indices(micro_edge, 1, dof_indices);
          extract_dof(dof_indices, flow, local_data);
          auto boundary_values_A = fe.evaluate_dof_at_boundary_points(local_data);

          filecsv << pressure(boundary_values_A.left) << "," << pressure(boundary_values_A.right);
        }
        filecsv << std::endl;
      }

      // concentration
      {
        std::fstream filecsv;
        filecsv.open(get_name("c", eid), std::ios::app);
        for (std::size_t micro_edge = 0; micro_edge < local_dof_map_transport.num_micro_edges(); micro_edge += 1) {
          if (micro_edge > 0)
            filecsv << ",";

          local_dof_map_transport.dof_indices(micro_edge, 0, dof_indices);
          extract_dof(dof_indices, transport, local_data);
          auto boundary_values = fe.evaluate_dof_at_boundary_points(local_data);

          filecsv << boundary_values.left << "," << boundary_values.right;
        }
        filecsv << std::endl;
      }
    }
  }
}

} // namespace macrocirculation