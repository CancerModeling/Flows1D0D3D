////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_CSV_WRITER_HPP
#define TUMORMODELS_GRAPH_CSV_WRITER_HPP

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
                 std::shared_ptr<GraphStorage> graph,
                 std::shared_ptr<DofMap> dof_map,
                 std::vector<std::string> component_names);

  void add_data(
    const std::shared_ptr< GraphStorage > & graph,
    const std::shared_ptr< DofMap > & dof_map,
    const std::vector< std::string >& component_names,
    const PetscVec& u);

  void add_data(
    const std::shared_ptr< GraphStorage > & graph,
    const std::shared_ptr< DofMap > & dof_map,
    const std::vector< std::string >& component_names,
    const std::vector< double >& u);

  void write(double t) const;

private:
  MPI_Comm d_comm;
  std::string d_foldername;
  std::string d_datasetname;
  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dof_map;
  std::vector<std::string> d_component_names;

  std::string get_name(std::size_t component, std::size_t vessel) const;
  std::string get_name_times() const;

  template < typename VectorType >
  struct Data {
    std::shared_ptr< GraphStorage > graph;
    std::shared_ptr< DofMap > dof_map;
    std::vector< std::string > component_names;
    const VectorType& u;
  };

  std::vector< Data< PetscVec > > d_petsc_data;
  std::vector< Data< std::vector< double > > > d_gmm_data;

  template <typename VectorType >
  void add_data_generic(
    const std::shared_ptr< GraphStorage > & graph,
    const std::shared_ptr< DofMap > & dof_map,
    const std::vector< std::string >& component_names,
    const VectorType& u);
};

} // namespace macrocirculation

#endif