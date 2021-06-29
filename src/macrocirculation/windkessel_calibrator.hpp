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
class ExplicitNonlinearFlowSolver;
class ImplicitLinearFlowSolver;

struct FlowData {
  /*! @brief Maps vertex ids to total flows. */
  std::map<size_t, double> flows;

  /*! @brief The sum of all flows. */
  double total_flow;
};

/*@! @brief Integrates the flow over the vessel tips for a single graph. */
class FlowIntegrator {
public:
  explicit FlowIntegrator(std::shared_ptr<GraphStorage> graph);

  void reset();

  /*! @brief Adds the flow contributions of the current time step to the total amount */
  void update_flow(const ExplicitNonlinearFlowSolver &solver, double tau);

  /*! @brief Adds the flow contributions of the current time step to the total amount */
  void update_flow(const ImplicitLinearFlowSolver &solver, double tau);

  /*! @brief Returns a data structure containing _all_ the flows and the total flow. */
  FlowData get_free_outflow_data() const;

protected:
  std::shared_ptr<GraphStorage> d_graph;

  /*! @brief Maps vertex ids to the total flow integrated over time. */
  std::map<size_t, double> d_total_flows;

  template<typename Solver>
  void update_flow_abstract(const Solver &solver, double tau);
};

double get_total_edge_capacitance(const std::shared_ptr<GraphStorage> &graph);

double get_total_edge_capacitance(const std::vector<std::shared_ptr<GraphStorage>> &list);

double get_total_flow(const std::vector<FlowData> &flows);

struct RCRData {
  double resistance;
  double capacitance;
};

class RCREstimator {
public:
  RCREstimator(std::vector<std::shared_ptr<GraphStorage>> graph_list, bool verbose);

  void reset();

  std::map<size_t, RCRData> estimate_parameters(const std::map<size_t, double> &total_flows, double total_flow);

  void set_total_C(double C);

  void set_total_R(double R);

protected:
  std::vector<std::shared_ptr<GraphStorage>> d_graph_list;

  /*! @brief Maps vertex ids to the total flow integrated over time. */
  std::map<size_t, double> d_total_flows;

  /*! @brief The total capacitance contributed by all the edges. */
  double d_total_C_edge;

  /*! @brief The total resistance of the whole human body. Per default in SI-Units. */
  double d_total_R;

  /*! @brief The total capacitance of the whole human body. Per default in SI-Units */
  double d_total_C;
};

void parameters_to_json(const std::string &filepath, const std::map<size_t, RCRData> &parameters, const std::shared_ptr<GraphStorage> &storage);

} // namespace macrocirculation


#endif //TUMORMODELS_WINDKESSEL_CALIBRATOR_HPP
