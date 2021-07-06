////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_CSV_WRITER_HPP
#define TUMORMODELS_GRAPH_CSV_WRITER_HPP

#include <map>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace macrocirculation {

class GraphStorage;
class DofMap;
class PetscVec;

class GraphCSVWriter {
public:
  GraphCSVWriter(MPI_Comm comm,
                 std::string foldername,
                 std::string datasetname,
                 const std::shared_ptr<GraphStorage> &graph);

  void add_setup_data(
    const std::shared_ptr<DofMap> &dof_map,
    size_t component_idx,
    const std::string &component_name);

  void setup();

  void add_data(const std::string &name, const PetscVec &u);
  void add_data(const std::string &name, const std::vector<double> &u);

  void write(double t);

private:
  MPI_Comm d_comm;
  std::string d_foldername;
  std::string d_datasetname;

  std::shared_ptr<GraphStorage> d_graph;

  bool is_setup;

  std::vector<double> time_steps;

  struct Data {
    std::shared_ptr<DofMap> dof_map;
    size_t component_idx;
    std::string component_name;
  };

  std::map<std::string, Data> d_data_map;

  std::map<std::string, std::reference_wrapper<const PetscVec>> petsc_data;
  std::map<std::string, std::reference_wrapper<const std::vector<double>>> gmm_data;

private:
  void reset_data();

  void clear_files();

  void write_meta_file();

  std::string get_meta_file_name() const;
  std::string get_csv_file_name(const std::string &component_name, size_t edge_id) const;

  template<typename VectorType>
  void write_generic(const Data &data, const VectorType &v) const;
};

} // namespace macrocirculation

#endif