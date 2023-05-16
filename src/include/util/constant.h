#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include <cmath>
#include <limits>

namespace xfdtd::constant {

constexpr double PI{3.14159265};
constexpr double EPS_0{8.854187817e-12};
constexpr double EPSILON_0{8.854187817e-12};
constexpr double MU_0{4 * PI * 1e-7};
constexpr double SQUARED_C_0{1 / (EPS_0 * MU_0)};
const double C_0{std::sqrt(SQUARED_C_0)};
constexpr double TOLERABLE_EPSILON{10 * std::numeric_limits<double>::epsilon()};
}  // namespace xfdtd::constant

#endif  // _CONSTANT_H_
