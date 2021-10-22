////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_RCR_ESTIMATOR_HPP
#define TUMORMODELS_RCR_ESTIMATOR_HPP

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

/*! @brief Integrates the flow over the vessel tips for a single graph. */
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

/*! @returns The total edge capacitance by summing over all ranks. */
double get_total_edge_capacitance(const std::shared_ptr<GraphStorage> &graph);

/*! @returns The total edge capacitance for a list of vessel networks by summing over all ranks. */
double get_total_edge_capacitance(const std::vector<std::shared_ptr<GraphStorage>> &list);

/*! @returns The total flow through all free-outflow boundaries. */
double get_total_flow(const std::vector<FlowData> &flows);

/*! @brief Simple data structure for the RCR system.
 *         Note that R1 and R2 are implicitly known, since R1 is the resistance of the 1D edge,
 *         while R2 = resistance - R1.
 */
struct RCRData {
  double resistance;
  double capacitance;
};

/*! @brief Estimates the quantities of the RCR-system for several coupled graphs from a list of total flows. */
class RCREstimator {
public:
  /*! @brief Constructs the RCR-estimator from a list of vessel networks. */
  RCREstimator(std::vector<std::shared_ptr<GraphStorage>> graph_list);

  /*! @brief Resets the total edge capacitance. */
  void reset();

  /*! @brief Estimates the RCR-data on the vessel tips for a single graph from the total flows at the vessel tips and the total flow through the whole network.
   *
   * @param total_flows A list of total flows through the free-outflow vessel tips of a _single_ graph.
   * @param total_flow The total flow through all free-outflow vessel tips for _all_ the graphs given in the graph_list of the constructor.
   * @return Maps vertex ids of the vessel tips to RCRData for the surrogate systems.
   */
  std::map<size_t, RCRData> estimate_parameters(const std::map<size_t, double> &total_flows, double total_flow);

  /*! @brief Sets the total (known) capacitance of the whole body/network. */
  void set_total_C(double C);

  /*! @brief Sets the total (known) resistance of the whole body/network. */
  void set_total_R(double R);

  double resistance_to_distribute() const;
  double capacitance_to_distribute() const;

  static double resistance_from_flow(double total_R, double flow);
  static double capacitance_from_flow(double total_C, double flow);

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


#endif //TUMORMODELS_RCR_ESTIMATOR_HPP
