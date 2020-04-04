////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_SYSTEMS_H
#define NETFV_SYSTEMS_H

// systems
#include "tumor.hpp"
#include "pressure.hpp"
#include "nutrient.hpp"
#include "hypoxic.hpp"
#include "necrotic.hpp"
#include "taf.hpp"
#include "grad_taf.hpp"
#include "ecm.hpp"
#include "mde.hpp"
#include "velocity.hpp"

namespace netfv {

/*!
 * @brief Compute matrix contribution from diffusion
 *
 * @param sys System for which assembly is performed
 */
void assemble_diffusion(BaseAssembly &sys);

/*!
 * @brief Compute matrix contribution from diffusion and advection
 *
 * @param sys System for which assembly is performed
 * @param pres Pressure system
 * @param tum Tumor system
 */
void assemble_diffusion_advection(BaseAssembly &sys, PressureAssembly &pres,
                                  TumAssembly &tum);

/*!
 * @brief Compute matrix contribution from advection
 *
 * @param sys System for which assembly is performed
 * @param pres Pressure system
 * @param tum Tumor system
 */
void assemble_advection(BaseAssembly &sys, PressureAssembly &pres, TumAssembly &tum);

} // namespace netfv

#endif // NETFV_SYSTEMS_H
