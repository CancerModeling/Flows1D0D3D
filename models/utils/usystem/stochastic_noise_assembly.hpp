////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#ifndef TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP
#define TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP

#include "usystem/abstraction.hpp"

namespace util {

/*! @brief Creates and assembles noise from a cylindrical Wiener process. */
class StochasticNoiseAssembly {
public:
  StochasticNoiseAssembly(unsigned int num_terms_bound,
                          unsigned int seed,
                          double scale,
                          double length,
                          double lower_bound,
                          double upper_bound);

  /*! @brief Assembles noise from a cylindrical Wiener process for the current time step.
   *
   *  @param assembly        The system to which we add the Wiener process.
   *  @param assembly_total  The total (tumor) system. If its values are not within an interval, the noise is deactivated.
   */
  void assemble(BaseAssembly &assembly, BaseAssembly &assembly_total) const;

  /*! @brief Calculates new stochastic coefficients for scaling our L2 basis.
   *         Should be called after every time step.
   *
   * @param dt  The current time step size, which determines the standard deviation of our stochastic increments.
   */
  void calculate_new_stochastic_coefficients(double dt);

  /*! @brief Calculates new stochastic coefficients for scaling our L2 basis.
   *         Should be called after every time step.
   *
   * Suppose \f$ S = H \Delta W \f$ is the stochastic term, \f$\Delta W = W_{n+1} - W_n\f$, and \f$ H\f$ is the heaviside function.
   * We compute
   * - \f$ s = \frac{1}{|\Omega|} \int_{\Omega} H \Delta W dx \f$
   * - \f$ h = \frac{1}{|\Omega|} \int_{\Omega} H dx \f$
   * - \f$ sh = s / h \f$ if \f$ h > 0 \f$ otherwise \f$ sh = 0 \f$. This term is calculated in variable d_avg
   *
   * Average corrected stochastic term is then given by \f$ S_{new} = H \Delta W - H sh = H (\Delta W - sh) \f$.
   *
   * We can check that \f$ \int_{\Omega} S_{new} = 0 \f$.
   *
   * @param dt  The current time step size, which determines the standard deviation of our stochastic increments.
   * @param assembly Assembly object to compute the average of stochastic term over the domain
   * @param total_assembly Assembly object to compute the average of stochastic term over the domain
   */
  void calculate_new_stochastic_coefficients(double dt, BaseAssembly &assembly, BaseAssembly &assembly_total);

  /*! @brief Evaluates the sum
   *         sum_{i,j,k}^J sqrt(8/L) * beta_{i,j,k} * cos(i*pi*x/L) * cos(j*pi*y/L) * cos(k*pi*z/L)
   *         where beta_{i,j,k} are our stochastic weights.
   *         This function can be used, if the assembly should take place directly in another assembly function
    *
    * @param p                  The point at which we want to evaluate the Wiener process.
    * @param field_value        The value of the field to which we add the noise.
    * @param total_field_value  The total value of the (whole tumor) field. The noise is deactivated if its value
    *                           is not within a certain interval.
    * @return The value of the Wiener process at the given quadrature point.
    */
  double eval_eigenfunctions_at_quadrature_point(const Point &p, double field_value, double total_field_value) const;

private:
  /*! @brief The an upper bound for the number of eigenfunctions. */
  unsigned int d_num_terms_bound;

  /*! @brief Random number generator. */
  std::default_random_engine d_generator;

  /*! @brief Length of one side of the cube. */
  double d_length;

  /*! @brief Scaling parameter for the stochastic disturbance. */
  double d_scale;

  /*! @brief Lower bound of the concentration to get the stochastic effect. */
  double d_lower_bound;

  /*! @brief Upper bound of the concentration to get the stochastic effect. */
  double d_upper_bound;

  /*! @brief Saves all the (current) stochastic weights of the 3D equation. */
  std::vector<double> d_stochastic_coefficients;

  /*! @brief Cached the last assembled contributions to the right hand side. Thus avoiding reassembly in a nonlinear loop. */
  mutable std::unique_ptr<NumericVector<Number>> d_cached_rhs;

  /*! @brief Flag indicating if the right hand side must be reassembled. */
  mutable bool d_reassemble;

  /*! @brief Average of stochastic function over domain */
  double d_avg;

private:
  /*! @brief Assembles noise from a cylindrical Wiener process for the current time step.
   *  @param assembly        The system to which we add the Wiener process.
   *  @param assembly_total  The total (tumor) system. If its values are not within an interval, the noise is deactivated.
   *  @param rhs             The right hand side vector to which we add the noise
   * */
  void assemble(BaseAssembly &assembly, BaseAssembly &assembly_total, NumericVector<Number> &rhs) const;

  std::size_t num_basis_functions() const;
};

} // namespace util

#endif //TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP
