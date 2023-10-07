#ifndef _FLOAT_COMPARE_H_
#define _FLOAT_COMPARE_H_

#include <cmath>

#include "util/constant.h"

namespace xfdtd {
inline bool isEqual(double a, double b,
                    double eps = constant::TOLERABLE_EPSILON) {
  return std::abs(a - b) <= eps;
}

inline bool isLessOrEqual(double a, double b,
                          double eps = constant::TOLERABLE_EPSILON) {
  return (std::abs(a - b) < eps) ? (true) : (a < b);
}

inline bool isGreaterOrEqual(double a, double b,
                             double eps = constant::TOLERABLE_EPSILON) {
  // fixed: Using integer absolute value function 'abs' when argument is of
  // floating point type
  return (std::abs(a - b) < eps) ? (true) : (a > b);
}

}  // namespace xfdtd

#endif  // _FLOAT_COMPARE_H_
