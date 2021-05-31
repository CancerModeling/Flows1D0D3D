////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_WINDKESSEL_CALIBRATOR_HPP
#define TUMORMODELS_WINDKESSEL_CALIBRATOR_HPP

#include <map>
#include <memory>
#include <mpi.h>

namespace macrocirculation {

class GraphStorage;
template<size_t degree>
class ExplicitNonlinearFlowSolver;

struct RCRData {
  double R_1;
  double R_2;
  double C;
};

class WindkesselCalibrator {
public:
  WindkesselCalibrator(std::shared_ptr<GraphStorage> graph, bool verbose);

  void reset();

  void reset_total_flows();

  void calculate_total_edge_capacitance();

  /*! @brief Adds the flow contributions of the current time step to the total amount */
  template<size_t degree>
  void update_flow(const ExplicitNonlinearFlowSolver<degree> &solver, double tau) {
    for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(MPI_COMM_WORLD))) {
      auto v = d_graph->get_vertex(v_id);
      if (v->is_leaf())
        d_total_flows[v_id] += solver.get_flow_at_vessel_tip(*v) * tau;
    }
  }

  std::map<size_t, RCRData> estimate_parameters();

  void set_total_C(double C);

  void set_total_R(double R);

protected:
  std::shared_ptr<GraphStorage> d_graph;

  bool d_verbose;

  /*! @brief Maps vertex ids to the total flow integrated over time. */
  std::map<size_t, double> d_total_flows;

  /*! @brief The total capacitance contributed by all the edges. */
  double d_total_C_edge;

  /*! @brief The total resistance of the whole human body. Per default in SI-Units. */
  double d_total_R;

  /*! @brief The total capacitance of the whole human body. Per default in SI-Units */
  double d_total_C;
};

} // namespace macrocirculation


#endif //TUMORMODELS_WINDKESSEL_CALIBRATOR_HPP
