////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFV_SYSTEMS_H
#define NETFV_SYSTEMS_H

// systems
#include "ecm.hpp"
#include "grad_taf.hpp"
#include "hypoxic.hpp"
#include "mde.hpp"
#include "necrotic.hpp"
#include "nutrient.hpp"
#include "pressure.hpp"
#include "taf.hpp"
#include "tumor.hpp"
#include "velocity.hpp"

namespace netfv {

// forward declare
class Model;

/*! @brief Compute matrix contribution from diffusion */
void assemble_diffusion(util::BaseAssembly &sys, Model *model);

/*! @brief Compute matrix contribution from diffusion and advection */
void assemble_diffusion_advection(util::BaseAssembly &sys, PressureAssembly &pres,
                                  TumAssembly &tum, Model *model);

/*! @brief Compute matrix contribution from advection */
void assemble_advection(util::BaseAssembly &sys, PressureAssembly &pres, TumAssembly &tum, Model *model);

} // namespace netfv

#endif // NETFV_SYSTEMS_H
