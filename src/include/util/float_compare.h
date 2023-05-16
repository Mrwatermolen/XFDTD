#ifndef _FLOAT_COMPARE_H_
#define _FLOAT_COMPARE_H_

#include <cassert>
#include <cmath>

namespace xfdtd {
inline bool isEqual(double a, double b, double eps) {
  return abs(a - b) <= eps;
}

inline bool isLessOrEqual(double a, double b, double eps) {
  return (abs(a - b) < eps) ? (true) : (a < b);
}

inline bool isGreaterOrEqual(double a, double b, double eps) {
  return (abs(a - b) < eps) ? (true) : (a > b);
}
}  // namespace xfdtd

#endif  // _FLOAT_COMPARE_H_