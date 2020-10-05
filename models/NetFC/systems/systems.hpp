////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef NETFC_SYSTEMS_H
#define NETFC_SYSTEMS_H

#include "../inp/inp.hpp"
#include "utilLibs.hpp"
#include "utils.hpp"

#include <string>

namespace netfc {

// forward declare Model class (needed for network system assembly)
class Model;

/*!
 * @brief Class to perform assembly
 *
 * Since we have two different equation system, one on 2d/3d mesh, and other
 * on 1d mesh, we can not use assemble functions which pass the equation
 * system to which they are attached.
 *
 * Here we store the pointer to Model object from which we can get both
 * equation systems easily.
 */
class ModelAssembly : public System::Assembly {

public:
  /*!
   * @brief Constructor
   */
  explicit ModelAssembly(Model *model) : d_model_p(model) {}

  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble() {}

private:
  /*! @brief Store the pointer to Model object */
  Model *d_model_p;
};

/**
 * @name Nutrient methods
 */
/**@{*/
class NutAssembly : public System::Assembly {
public:
  NutAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_nut(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

void boundary_condition_nut(EquationSystems &es);

/** @}*/

/**
 * @name TAF methods
 */
/**@{*/
class TafAssembly : public System::Assembly {
public:
  TafAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_taf(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/** @}*/

/**
 * @name Grad TAF methods
 */
/**@{*/
class GradTafAssembly : public System::Assembly {
public:
  GradTafAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_2d();
  void assemble_3d();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_grad_taf(const Point &p, const Parameters &es,
                                  const std::string &system_name,
                                  const std::string &var_name);

/** @}*/

/**
 * @name Hypoxic methods
 */
/**@{*/
class HypAssembly : public System::Assembly {
public:
  HypAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_hyp(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

Number initial_condition_hyp_kernel(const Point &p,
                                    const netfc::InputDeck *deck);

/** @}*/

/**
 * @name Necrotic methods
 */
/**@{*/
class NecAssembly : public System::Assembly {
public:
  NecAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_nec(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/** @}*/

/**
 * @name Tumor methods
 */
/**@{*/
class TumAssembly : public System::Assembly {
public:
  TumAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_tum(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

/** @}*/

/**
 * @name ECM methods
 */
/**@{*/
class EcmAssembly : public System::Assembly {
public:
  EcmAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_ecm(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

void boundary_condition_ecm(EquationSystems &es);

/** @}*/

/**
 * @name MDE methods
 */
/**@{*/
class MdeAssembly : public System::Assembly {
public:
  MdeAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  void assemble_2();
  void assemble_3();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_mde(const Point &p, const Parameters &es,
                             const std::string &system_name, const std::string &var_name);

void boundary_condition_mde(EquationSystems &es);

/** @}*/

/**
 * @name Pressure
 */
/**@{*/

class PressureAssembly : public System::Assembly {
public:
  PressureAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  void assemble_1();
  std::string d_sys_name;
  Model *d_model_p;
};

Number initial_condition_pres(const Point &p, const Parameters &es,
                              const std::string &system_name, const std::string &var_name);

void boundary_condition_pres(EquationSystems &es);

class VelocityAssembly : public System::Assembly {
public:
  VelocityAssembly(Model *model, const std::string system_name) : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;
  std::string d_sys_name;
  Model *d_model_p;

private:
  void assemble_2d();
  void assemble_3d();
};

/** @}*/

} // namespace netfc

#endif // NETFC_SYSTEMS_H
