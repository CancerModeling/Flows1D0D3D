#ifndef TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP
#define TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP

#include "usystem/abstraction.hpp"

namespace util {

/*! @brief Creates and assembles noise from a cylindrical Wiener process. */
class StochasticNoiseAssembly {
public:
  StochasticNoiseAssembly(unsigned int num_eigenfunctions,
                          unsigned int seed,
                          double scale,
                          double length,
                          double nutrient_lower_bound,
                          double nutrient_upper_bound);

  /*! @brief Assembles noise from a cylindrical Wiener process for the current timestep. */
  void assemble(BaseAssembly &assembly) const;

  /*! @brief Calculates new stochastic coefficients for scaling our L2 basis.
   *         Should be called after every timestep.
   *
   * @param dt  The current time step size, which determines the standard deviation of our stochastic increments.
   */
  void calculate_new_stochastic_coefficients(double dt);

private:
  /*! @brief The number of eigenfuntions of our covariance matrix. */
  unsigned int d_num_eigenfunctions;

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

private:
  /*! @brief Evaluates the sum
   *         sum_{i,j,k}^J sqrt(8/L) * beta_{i,j,k} * cos(i*pi*x/L) * cos(j*pi*y/L) * cos(k*pi*z/L)
   *         where beta_{i,j,k} are our stochastic weights.
   */
  double eval_eigenfunctions_at_quadrature_point(const Point &p, double field_value) const;
};

} // namespace util

#endif //TUMORMODELS_STOCHASTIC_NOISE_ASSEMBLY_HPP
