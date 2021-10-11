////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_VESSEL_TREE_FLOW_INTEGRATOR_HPP
#define TUMORMODELS_VESSEL_TREE_FLOW_INTEGRATOR_HPP

#include "mpi.h"
#include "graph_storage.hpp"

namespace macrocirculation {

class PetscVec;
class DofMap;

struct VesselTreeFlowIntegratorResult {
  /*! @brief Coordinates of the vessel tip. */
  Point point;

  /*! @brief Vertex id at the vessel tip. */
  size_t vertex_id;

  /*! @brief Number of levels at the vessel tip. */
  size_t levels;

  /*! @brief the averaged flow value in [cm^3/s]. */
  double averaged_flow;

  /*! @brief the averaged 3D pressure value in [Ba]. */
  double averaged_3D_pressure;
};

class VesselTreeFlowIntegrator {
public:
  VesselTreeFlowIntegrator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap> dof_map);

  void reset();

  void add(const PetscVec &u, double tau);

  std::vector<VesselTreeFlowIntegratorResult> calculate();

  void write(const std::string &folder, const std::string &dataset_name);

  static void write(const std::string &folder, const std::string &dataset_name, const GraphStorage& graph, const std::vector<VesselTreeFlowIntegratorResult> &results);

private:
  MPI_Comm d_comm;
  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dof_map;

  struct AvgVesselTipData {
    double flow;
    double pressure_1d;
    double pressure_3d;
    double time;
  };

  std::map<size_t, AvgVesselTipData> d_avg_data;
};

}

#endif //TUMORMODELS_VESSEL_TREE_FLOW_INTEGRATOR_HPP
