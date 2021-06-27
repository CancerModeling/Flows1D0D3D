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

struct RCRData {
  double resistance;
  double capacitance;
};

class WindkesselCalibrator {
public:
  WindkesselCalibrator(std::shared_ptr<GraphStorage> graph, bool verbose);

  void reset();

  void reset_total_flows();

  void calculate_total_edge_capacitance();

  /*! @brief Adds the flow contributions of the current time step to the total amount */
  void update_flow(const ExplicitNonlinearFlowSolver &solver, double tau);

  std::map<size_t, RCRData> estimate_parameters_local();

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

void parameters_to_json(const std::string &filepath, const std::map<size_t, RCRData> &parameters, const std::shared_ptr<GraphStorage> &storage);

} // namespace macrocirculation


#endif //TUMORMODELS_WINDKESSEL_CALIBRATOR_HPP
