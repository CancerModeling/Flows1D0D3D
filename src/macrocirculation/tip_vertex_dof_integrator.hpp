////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_TIP_VERTEX_DOF_INTEGRATOR_HPP
#define TUMORMODELS_TIP_VERTEX_DOF_INTEGRATOR_HPP

#include <memory>
#include <vector>
#include <map>

#include "mpi.h"

namespace macrocirculation {

// forward declarations
class GraphStorage;
class ExplicitNonlinearFlowSolver;
class ImplicitLinearFlowSolver;
class DofMap;
class PetscVec;

/*! @brief Integrates the vertex dofs over time. */
class TipVertexDofIntegrator {
public:
  explicit TipVertexDofIntegrator(MPI_Comm comm, std::shared_ptr<GraphStorage> graph, std::shared_ptr<DofMap>);

  void reset();

  /*! @brief Adds the flow contributions of the current time step to the total amount */
  void update_vertex_dof(const PetscVec &vec, double tau);

  /*! @brief Time interval over which the integrator was applied. */
  double get_integration_time() const;

  /*! @brief Returns a data structure containing the integrated vertex dof values specified in the given vector. */
  std::map<size_t, std::vector< double > > get_integral_value(const std::vector< size_t >& vertex_dof_numbers) const;

protected:
  MPI_Comm d_comm;

  std::shared_ptr<GraphStorage> d_graph;

  std::shared_ptr<DofMap> d_dof_map;

  double d_total_integration_time;

  /*! @brief Map vertex ids to vectors containing the average quantity. */
  std::map< size_t, std::vector< double > > d_quantities;
};

}

#endif //TUMORMODELS_TIP_VERTEX_DOF_INTEGRATOR_HPP