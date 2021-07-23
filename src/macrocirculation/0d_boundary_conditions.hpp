////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H
#define TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H

#include <memory>
#include <string>

namespace macrocirculation {

// forward declarations:
class GraphStorage;

/*! @brief Sets the boundary conditions using a tree estimate. */
void set_0d_tree_boundary_conditions(const std::shared_ptr<GraphStorage> &graph, const std::string &name_prefix);

/*! @brief Takes the RCR boundary conditions and partitions the conductances and resistances into an RCR system.
 *         Every pressure corresponds to the mean pressure of a part of the vascular tree.
 *         Currently the partitioning is as follows:
 *         - 4 RCR-systems for the arterioles,
 *         - 1 RCR-system for the capillaries,
 *         - 1 RCR-system for the venlues and
 *         - 1 RCR-system for the small veins
 *
 * @param graph The vessel network whose boundary conditions are changed.
 */
void convert_rcr_to_partitioned_tree_bcs(const std::shared_ptr<GraphStorage> &graph);

void convert_rcr_to_rcl_chain_bcs(const std::shared_ptr<GraphStorage> &graph);

} // namespace macrocirculation

#endif //TUMORMODELS_SET_0D_TREE_BOUNDARY_CONDITIONS_H
