////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_GRAPH_FLOW_AND_CONCENTRATION_WRITER_HPP
#define TUMORMODELS_GRAPH_FLOW_AND_CONCENTRATION_WRITER_HPP

#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace macrocirculation {

class GraphStorage;
class DofMap;
class PetscVec;

class GraphFlowAndConcentrationWriter {
public:
  GraphFlowAndConcentrationWriter(MPI_Comm comm,
                                  std::string foldername,
                                  std::string datasetname,
                                  std::shared_ptr<GraphStorage> graph,
                                  std::shared_ptr<DofMap> dof_map_flow,
                                  std::shared_ptr<DofMap> dof_map_transport);

  void write(double t, const std::vector<double> &flow, const std::vector<double> &transport) const;

  void write(double t, const PetscVec &flow, const PetscVec &transport) const;

private:
  MPI_Comm d_comm;
  std::string d_foldername;
  std::string d_datasetname;
  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dof_map_flow;
  std::shared_ptr<DofMap> d_dof_map_transport;

  std::string get_name(const std::string &name_quantity, std::size_t vessel) const;
  std::string get_name_times() const;

  template < typename GenericVectorType >
  void write_generic(double t, const GenericVectorType &flow, const GenericVectorType &transport) const;
};

} // namespace macrocirculation

#endif