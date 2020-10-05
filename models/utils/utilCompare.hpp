////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_COMPARE_H
#define UTIL_COMPARE_H

#include <cmath>

// Assign tolerance for comparison
#define COMPARE_EPS 1e-5

namespace util {

/*!
 * @brief Compares if a is approximately equal to b
 * @param a Value a
 * @param b Value b
 * @return Result true if approximately equal else false
 */
inline bool approximatelyEqual(const double &a, const double &b) {
  return std::abs(a - b) <=
         ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
          COMPARE_EPS);
}

/*!
 * @brief Compares if a is essentially equal to b
 * @param a Value a
 * @param b Value b
 * @return Result true if essentially equal else false
 */
inline bool essentiallyEqual(const double &a, const double &b) {
  return std::abs(a - b) <=
         ((std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a)) *
          COMPARE_EPS);
}

/*!
 * @brief Compares if a > to b
 * @param a Value a
 * @param b Value b
 * @return Result true if a is definitely greater than b
 */
inline bool definitelyGreaterThan(const double &a, const double &b) {
  return (a - b) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
                    COMPARE_EPS);
}


/*!
 * @brief Compares if a is < to b
 * @param a Value a
 * @param b Value b
 * @return Result true if a is definitely less than b
 */
inline bool definitelyLessThan(const double &a, const double &b) {
  return (b - a) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
                    COMPARE_EPS);
}

} // namespace util

#endif // UTIL_COMPARE_H