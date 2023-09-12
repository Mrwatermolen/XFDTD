#ifndef _TYPE_DEFINE_H_
#define _TYPE_DEFINE_H_

#include <algorithm>
#include <memory>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>

namespace xfdtd {
class Source;
class Object;
class Boundary;
class YeeCell;
class Monitor;
// coordinate
using PointVector = xt::xtensor_fixed<double, xt::xshape<3>>;

// using DoubleArrary1D = std::vector<double>;
// using DoubleArrary2D = Eigen::MatrixXd;
// using DoubleArrary3D = Eigen::Tensor<double, 3>;
using DoubleArrary1D = xt::xarray<double>;
using DoubleArrary2D = xt::xarray<double>;
using DoubleArrary3D = xt::xtensor<double, 3>;

// using ElectromagneticFieldTimedomainArray = Eigen::Tensor<double, 3>;
// using ElectromagneticFieldFrequencydomainArray =
//     Eigen::Tensor<std::complex<double>, 3>;
using ElectromagneticFieldTimedomainArray = xt::xarray<double>;
using ElectromagneticFieldFrequencydomainArray =
    xt::xarray<std::complex<double>>;
using EFTA = ElectromagneticFieldTimedomainArray;
using EFFA = ElectromagneticFieldFrequencydomainArray;

using SpatialIndex = int;

using ObjectArray = std::vector<std::shared_ptr<Object>>;
using BoundaryArray = std::vector<std::shared_ptr<Boundary>>;
using YeeCellArray = std::vector<std::shared_ptr<YeeCell>>;
using MonitorArray = std::vector<std::shared_ptr<Monitor>>;

enum class PlaneType { XY, YZ, ZX };
enum class CoordinateNormalVector { XN, YN, ZN, XP, YP, ZP };
using CNV = CoordinateNormalVector;

inline DoubleArrary3D allocateDoubleArray3D(SpatialIndex nx, SpatialIndex ny,
                                            SpatialIndex nz,
                                            double default_value = 0.0) {
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

inline DoubleArrary2D allocateDoubleArray2D(SpatialIndex nx, SpatialIndex ny,
                                            double default_value = 0.0) {
  std::vector<size_t> a_shape = {static_cast<size_t>(nx),
                                 static_cast<size_t>(ny)};
  auto arr{DoubleArrary2D(a_shape)};
  arr.fill(default_value);
  return arr;
}

inline DoubleArrary1D allocateDoubleArray1D(SpatialIndex n,
                                            double default_value = 0.0) {
  std::vector<size_t> a_shape = {static_cast<size_t>(n)};
  auto arr{DoubleArrary1D(a_shape)};
  arr.fill(default_value);
  return arr;
}
}  // namespace xfdtd

#endif  // _TYPE_DEFINE_H_
