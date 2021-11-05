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
#include <nlohmann/json.hpp>
#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"

namespace macrocirculation {

GraphCSVWriter::GraphCSVWriter(MPI_Comm comm,
                               std::string foldername,
                               std::string datasetname,
                               std::shared_ptr<GraphStorage> graph)
    : d_comm(comm),
      d_foldername(std::move(foldername)),
      d_datasetname(std::move(datasetname)),
      d_graph(std::move(graph)),
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
  write_meta_file();
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
  write_times();

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

void GraphCSVWriter::write_times() {
  if (mpi::rank(d_comm) == 0) {
    std::ofstream f(get_time_csv_file_path(), std::ios::out);

    for (size_t k = 0; k < time_steps.size(); k += 1) {
      if (k != 0)
        f << ", ";
      f << time_steps[k];
    }
  }
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
      std::ofstream(get_csv_file_path(name, eid), std::ios::out);
  }

  if (mpi::rank(d_comm) == 0) {
    std::ofstream(get_time_csv_file_path(), std::ios::out);
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
      std::fstream filecsv(get_csv_file_path(data.component_name, eid), std::ios::app);

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
  name << d_datasetname << "_" << component_name << "_vessel" << std::setfill('0') << std::setw(5) << edge_id << ".csv";
  return name.str();
}

std::string GraphCSVWriter::get_csv_file_path(const std::string &component_name, size_t edge_id) const {
  return d_foldername + "/" + get_csv_file_name(component_name, edge_id);
}

std::string GraphCSVWriter::get_time_csv_file_name() const {
  return d_datasetname + "_time.csv";
}

std::string GraphCSVWriter::get_time_csv_file_path() const {
  return d_foldername + "/" + get_time_csv_file_name();
}

void GraphCSVWriter::write_meta_file() {
  // only one rank writes the file
  if (mpi::rank(d_comm) != 0)
    return;

  std::ofstream f(get_meta_file_name(), std::ios::out);

  using json = nlohmann::json;

  json j;

  auto vessel_list = json::array();
  for (auto eid : d_graph->get_edge_ids()) {
    const auto edge = d_graph->get_edge(eid);
    const double length = edge->has_physical_data() ? edge->get_physical_data().length : 1.;
    const auto h = length / static_cast<double>(edge->num_micro_edges());

    const auto &vertex_left = *d_graph->get_vertex(edge->get_vertex_neighbors()[0]);
    const auto &vertex_right = *d_graph->get_vertex(edge->get_vertex_neighbors()[1]);

    std::vector<double> coordinates;
    for (std::size_t micro_edge = 0; micro_edge < edge->num_micro_edges(); micro_edge += 1) {
      coordinates.push_back(h * micro_edge);
      coordinates.push_back(h * (micro_edge + 1));
    }

    auto &pdata = edge->get_physical_data();

    json filepath_obj;
    for (auto d : d_data_map) {
      filepath_obj[d.first] = get_csv_file_name(d.first, eid);
    }

    json vertices_obj = {
      {"left",
       {{"id", vertex_left.get_id()},
        {"name", vertex_left.get_name()}}},
      {"right",
       {{"id", vertex_right.get_id()},
        {"name", vertex_right.get_name()}}}};

    json vessel_obj = {
      {"edge_id", edge->get_id()},
      {"name", edge->get_name()},
      {"coordinates", coordinates},
      {"filepaths", filepath_obj},
      {"vertices", vertices_obj},
      {"A0", pdata.A0},
      {"G0", pdata.G0}};

    vessel_list.push_back(vessel_obj);
  }

  j["vessels"] = vessel_list;
  j["filepath_time"] = get_time_csv_file_name();

  f << j.dump(1);
}

} // namespace macrocirculation