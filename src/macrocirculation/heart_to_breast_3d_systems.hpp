////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_HEART_TO_BREAST_3_D_SYSTEMS_HPP
#define TUMORMODELS_HEART_TO_BREAST_3_D_SYSTEMS_HPP

#include "libmesh/libmesh.h"
#include "assembly_system.hpp"

#include <memory>


namespace macrocirculation {

// forward declarations:
class HeartToBreast3DSolver;

/*! @brief Initial condition function for pressure systems. */
inline lm::Number ic_p(const lm::Point &p, const lm::Parameters &es,
                const std::string &system_name,
                const std::string &var_name) { return 40000.; }

/*! @brief Initial condition function for tumor system. */
inline lm::Number ic_tum(const lm::Point &p, const lm::Parameters &es,
                       const std::string &system_name,
                       const std::string &var_name) { return 0.; }

/*! @brief Function that applies initial condition to pressure systems. */
inline void ic(lm::EquationSystems &es, const std::string &system_name) {
  auto &sys = es.get_system<lm::TransientLinearImplicitSystem>(
    system_name);
  if (system_name == "Capillary_Pressure" or system_name == "Tissue_Pressure")
    sys.project_solution(ic_p, nullptr, es.parameters);
}

/*! @brief Class that handles assembly of matrix and right-hand side of capillary pressure system. */
class CapillaryPressure : public BaseAssembly {
public:
  CapillaryPressure(HeartToBreast3DSolver *model, lm::MeshBase &mesh,
        lm::TransientLinearImplicitSystem &sys)
    : BaseAssembly("Capillary_Pressure", mesh, sys, 1,
                       {sys.variable_number("p_cap")}),
      d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  /*! @brief Assemble matrix and right-hand side. */
  void assemble() override;

  /*! @brief Assemble contribution to matrix and right-hand side from 1D system. */
  void assemble_1d();

  /*! @brief Pointer to the 3D model class to access relevant data and methods. */
  HeartToBreast3DSolver *d_model_p;
};

/*! @brief Class that handles assembly of matrix and right-hand side of tissue pressure system. */
class TissuePressure : public BaseAssembly {
public:
  TissuePressure(HeartToBreast3DSolver *model, lm::MeshBase &mesh,
                    lm::TransientLinearImplicitSystem &sys)
    : BaseAssembly("Tissue_Pressure", mesh, sys, 1,
                       {sys.variable_number("p_tis")}),
      d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  /*! @brief Assemble matrix and right-hand side. */
  void assemble() override;

  /*! @brief Pointer to the 3D model class to access relevant data and methods. */
  HeartToBreast3DSolver *d_model_p;
};

/*! @brief Class that handles assembly of matrix and right-hand side of capillary pressure system. */
class CapillaryNutrient : public BaseAssembly {
public:
  CapillaryNutrient(HeartToBreast3DSolver *model, lm::MeshBase &mesh,
                    lm::TransientLinearImplicitSystem &sys)
    : BaseAssembly("Capillary_Nutrient", mesh, sys, 1,
                   {sys.variable_number("nut_cap")}),
      d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  /*! @brief Assemble matrix and right-hand side. */
  void assemble() override;

  /*! @brief Assemble contribution to matrix and right-hand side from 1D system. */
  void assemble_1d();

  /*! @brief Pointer to the 3D model class to access relevant data and methods. */
  HeartToBreast3DSolver *d_model_p;
};

/*! @brief Class that handles assembly of matrix and right-hand side of tissue pressure system. */
class TissueNutrient : public BaseAssembly {
public:
  TissueNutrient(HeartToBreast3DSolver *model, lm::MeshBase &mesh,
                 lm::TransientLinearImplicitSystem &sys)
    : BaseAssembly("Tissue_Nutrient", mesh, sys, 1,
                   {sys.variable_number("nut_tis")}),
      d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  /*! @brief Assemble matrix and right-hand side. */
  void assemble() override;

  /*! @brief Pointer to the 3D model class to access relevant data and methods. */
  HeartToBreast3DSolver *d_model_p;
};

/*! @brief Class that handles assembly of matrix and right-hand side of tissue pressure system. */
class Tumor : public BaseAssembly {
public:
  Tumor(HeartToBreast3DSolver *model, lm::MeshBase &mesh,
                 lm::TransientLinearImplicitSystem &sys)
    : BaseAssembly("Tumor", mesh, sys, 2,
                   {sys.variable_number("tum"), sys.variable_number("mu_tum")}),
      d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  /*! @brief Assemble matrix and right-hand side. */
  void assemble() override;

  /*! @brief Pointer to the 3D model class to access relevant data and methods. */
  HeartToBreast3DSolver *d_model_p;
};

} // namespace macrocirculation


#endif //TUMORMODELS_HEART_TO_BREAST_3_D_SYSTEMS_HPP
