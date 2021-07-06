////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "graph_csv_writer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

GraphCSVWriter::GraphCSVWriter(MPI_Comm comm,
                               std::string foldername,
                               std::string datasetname,
                               const std::shared_ptr<GraphStorage> &graph)
    : d_comm(comm),
      d_foldername(std::move(foldername)),
      d_datasetname(std::move(datasetname)),
      d_graph(graph),
      is_setup(false) {}

void GraphCSVWriter::add_setup_data(
  const std::shared_ptr<DofMap> &dof_map,
  size_t component_idx,
  const std::string &component_name) {
  if (is_setup)
    throw std::runtime_error("cannot add setup data after setup was called");

  d_data_map[component_name] = {dof_map, component_idx, component_name};
}

void GraphCSVWriter::setup() {
  is_setup = true;
  clear_files();
}

void GraphCSVWriter::add_data(const std::string &name, const PetscVec &u) {
  if (!is_setup)
    throw std::runtime_error("data can only be added after calling setup");
  if (d_data_map.count(name) == 0)
    throw std::runtime_error("key " + name + " was not added in setup phase");
  petsc_data.emplace(name, std::cref(u));
}

void GraphCSVWriter::add_data(const std::string &name, const std::vector<double> &u) {
  if (!is_setup)
    throw std::runtime_error("data can only be added after calling setup");
  if (d_data_map.count(name) == 0)
    throw std::runtime_error("key " + name + " was not added in setup phase");
  gmm_data.emplace(name, std::cref(u));
}

void GraphCSVWriter::write(double t) {
  // add time step
  time_steps.push_back(t);
  write_meta_file();

  // write the added data to disk
  for (auto &data_it : d_data_map) {
    auto name = data_it.first;

    if (petsc_data.count(name)) {
      write_generic(data_it.second, petsc_data.at(name));
    } else if (gmm_data.count(name)) {
      write_generic(data_it.second, gmm_data.at(name));
    } else {
      throw std::runtime_error("component " + name + " was not added for writing");
    }
  }

  //
  reset_data();
}

void GraphCSVWriter::reset_data() {
  petsc_data.clear();
  gmm_data.clear();
}

void GraphCSVWriter::clear_files() {
  // all csv files get cleared
  for (auto &data_it : d_data_map) {
    auto name = data_it.first;

    for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm)))
      std::ofstream(get_csv_file_name(name, eid), std::ios::out);
  }
}

template<typename VectorType>
void GraphCSVWriter::write_generic(const Data &data, const VectorType &v) const {
  auto dof_map = data.dof_map;
  auto component = data.component_idx;

  // write vessel data
  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto &edge = d_graph->get_edge(eid);
    const auto &local_dof_map = dof_map->get_local_dof_map(*edge);
    const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
    const auto h = length / local_dof_map.num_micro_edges();

    FETypeNetwork fe(create_midpoint_rule(), local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions(), 0);
    std::vector<double> local_data(local_dof_map.num_basis_functions(), 0);

    // append line to csv file
    {
      std::fstream filecsv(get_csv_file_name(data.component_name, eid), std::ios::app);

      for (std::size_t micro_edge = 0; micro_edge < local_dof_map.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";

        local_dof_map.dof_indices(micro_edge, component, dof_indices);
        extract_dof(dof_indices, v, local_data);
        auto boundary_values = fe.evaluate_dof_at_boundary_points(local_data);

        filecsv << boundary_values.left << "," << boundary_values.right;
      }
      filecsv << std::endl;
    }
  }
}

std::string GraphCSVWriter::get_meta_file_name() const {
  return d_foldername + "/" + d_datasetname + ".json";
}

std::string GraphCSVWriter::get_csv_file_name(const std::string &component_name, size_t edge_id) const {
  std::stringstream name;
  name << d_foldername << "/" << d_datasetname << "_" << component_name << "_vessel" << std::setfill('0') << std::setw(5) << edge_id
       << ".csv";
  return name.str();
}

void GraphCSVWriter::write_meta_file() {
  // only one rank writes the file
  if (mpi::rank(d_comm) != 0)
    return;

  for (auto eid : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto edge = d_graph->get_edge(eid);
    const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
    const auto h = length / static_cast<double>(edge->num_micro_edges());

    std::vector<double> coordinates;
    for (std::size_t micro_edge = 0; micro_edge < edge->num_micro_edges(); micro_edge += 1) {
      coordinates.push_back(h * micro_edge);
      coordinates.push_back(h * (micro_edge + 1));
    }
  }
}

/*
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

    {
      std::fstream filecsv;
      std::stringstream name;
      name << d_foldername << "/" << d_datasetname << "_"
           << "p"
           << "_vessel" << std::setfill('0') << std::setw(5) << eid << ".csv";
      filecsv.open(name.str(), std::ios::out);
      for (std::size_t micro_edge = 0; micro_edge < local_dof_map.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";
        filecsv << h * micro_edge << "," << h * (micro_edge + 1);
      }
      filecsv << std::endl;
    }
  }

  if (mpi::rank(d_comm) == 0) {
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

void GraphCSVWriter::write(double t, const std::vector<double> &data) const {
  write_generic(t, data);
}

void GraphCSVWriter::write(double t, const PetscVec &data) const {
  write_generic(t, data);
}

template<typename FunctionType>
void GraphCSVWriter::write_generic(double time, const FunctionType &data) const {
  // write time step information
  if (mpi::rank(d_comm) == 0) {
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

    // TODO: make the pressure calculation clearer
    if (edge->has_physical_data() && d_component_names[0] == "Q" && d_component_names[0] == "A") {
      const auto A0 = edge->get_physical_data().A0;
      const auto G0 = edge->get_physical_data().G0;

      const auto pressure = [=](auto A) { return G0 * (std::sqrt(A / A0) - 1.0) / 1.33332; };

      std::fstream filecsv;
      std::stringstream name;
      name << d_foldername << "/" << d_datasetname << "_"
           << "p"
           << "_vessel" << std::setfill('0') << std::setw(5) << eid << ".csv";
      filecsv.open(name.str(), std::ios::app);
      for (std::size_t micro_edge = 0; micro_edge < local_dof_map.num_micro_edges(); micro_edge += 1) {
        if (micro_edge > 0)
          filecsv << ",";

        local_dof_map.dof_indices(micro_edge, 1, dof_indices);
        extract_dof(dof_indices, data, local_data);
        auto boundary_values_A = fe.evaluate_dof_at_boundary_points(local_data);

        filecsv << pressure(boundary_values_A.left) << "," << pressure(boundary_values_A.right);
      }
      filecsv << std::endl;
    }
  }
}
 */

} // namespace macrocirculation