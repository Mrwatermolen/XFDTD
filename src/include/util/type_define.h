#ifndef _TYPE_DEFINE_H_
#define _TYPE_DEFINE_H_

#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>

namespace xfdtd {
// coordinate
using PointVector = xt::xtensor_fixed<double, xt::xshape<3>>;

enum class Orientation { XN, XP, YN, YP, ZN, ZP };
enum class Axis { X, Y, Z };

class SizeInf {
 public:
  SizeInf() = default;
  ~SizeInf() = default;
};

using ElectromagneticFieldTimeDomainArray = xt::xarray<double>;
using ElectromagneticFieldFrequencyDomainArray =
    xt::xarray<std::complex<double>>;
using EFTA = ElectromagneticFieldTimeDomainArray;
using EFFA = ElectromagneticFieldFrequencyDomainArray;

using SpatialIndex = size_t;

enum class PlaneType { XY, YZ, ZX };

inline auto allocateDoubleArray3D(SpatialIndex nx, SpatialIndex ny,
                                  SpatialIndex nz, double default_value = 0.0) {
  std::vector<size_t> a_shape = {static_cast<size_t>(nx),
                                 static_cast<size_t>(ny),
                                 static_cast<size_t>(nz)};
  auto arr{
      xt::xtensor<double, 3>({static_cast<size_t>(nx), static_cast<size_t>(ny),
                              static_cast<size_t>(nz)},
                             default_value)};
  // arr.fill(default_value);
  return arr;
}

inline xt::xarray<double> allocateDoubleArray1D(SpatialIndex n,
                                                double default_value = 0.0) {
  std::vector<size_t> a_shape = {static_cast<size_t>(n)};
  auto arr{xt::xarray<double>(a_shape)};
  arr.fill(default_value);
  return arr;
}
}  // namespace xfdtd

#endif  // _TYPE_DEFINE_H_
