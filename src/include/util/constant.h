#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#include <cmath>
#include <limits>

namespace xfdtd::constant {

//
inline constexpr double PI{M_PI};

//
inline constexpr double EPSILON_0{8.854187817e-12};

//
inline constexpr double MU_0{4 * PI * 1e-7};

//
inline constexpr double SQUARED_C_0{1 / (EPSILON_0 * MU_0)};

// the intrinsic impedance of free space
inline const double ETA_0{std::sqrt(constant::MU_0 / constant::EPSILON_0)};

//
inline const double C_0{std::sqrt(SQUARED_C_0)};

//
inline constexpr double TOLERABLE_EPSILON{
    10 * std::numeric_limits<double>::epsilon()};
}  // namespace xfdtd::constant

#endif  // _CONSTANT_H_
