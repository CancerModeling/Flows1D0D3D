////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef PARTIONING_PERFUSION_HPP
#define PARTIONING_PERFUSION_HPP

#include <vector>
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "logger.hpp"
#include "libmesh_includes.hpp"

namespace lm = libMesh;

namespace macrocirculation {

/*! @brief Compute estimated flow rate
 * @param radii Vector of outlet radius
 * @param gamma Exponent to estimate flow rate
 * @param flow Vector of estimated flow rates for each outlet
 */
void estimate_flow_rate(const std::vector<double> &radii,
                                double gamma,
                                std::vector<double> &flow);

/*! @brief Compute estimated flow rate
 * @param radii Vector of outlet radius
 * @param flow Vector of estimated flow rates for each outlet
 * @param vol Vector of estimated volume fraction
 */
void estimate_vol_fraction(const std::vector<double> &radii,
                        const std::vector<double> &flow,
                        std::vector<double> &vol);

/*! @brief Creates piecewise constant field where value of field at each element is given by the partition data
 * @param mesh Libmesh mesh object
 * @param field Libmesh field
 * @param partition Vector of integers for each element
 */
void create_partition_field(const lm::MeshBase &mesh,
                            lm::ExplicitSystem &field,
                            const std::vector<long> &partition);

/*! @brief Given outlet points and radii of vessel at outlet points, this function partitions the mesh into different territories
 * @param x Vector of coordinates of outlet points
 * @param r Vector of outlet radius
 * @param v Volume fractions
 * @param mesh Libmesh mesh
 * @param partition Data that contains the outlet id for each element in the mesh
 * @param n_rand_loop Number of loop for redistribution of elements
 * @param log Logger
 */
void create_perfusion_territory(const std::vector<lm::Point> &x,
                                const std::vector<double> &v,
                                lm::ReplicatedMesh &mesh,
                                std::vector<long> &partition,
                                int n_rand_loop,
                                lm::ExplicitSystem &partition_field,
                                lm::EquationSystems &eq_sys,
                                macrocirculation::Logger &log,
                                std::string out_dir,
                                int out_freq);

} // namespace macrocirculation

#endif //PARTIONING_PERFUSION_HPP