#ifndef _TYPE_DEFINE_H_
#define _TYPE_DEFINE_H_

#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace xfdtd {
class Source;
class Object;
class Boundary;
class YeeCell;
class Monitor;
// coordinate
using PointVector = Eigen::Vector3d;

using DoubleArrary1D = std::vector<double>;
using DoubleArrary2D = Eigen::MatrixXd;
using DoubleArrary3D = Eigen::Tensor<double, 3>;

using ElectromagneticFieldTimedomainArray = Eigen::Tensor<double, 3>;
using ElectromagneticFieldFrequencydomainArray =
    Eigen::Tensor<std::complex<double>, 3>;
using EFTA = ElectromagneticFieldTimedomainArray;
using EFFA = ElectromagneticFieldFrequencydomainArray;

using SpatialIndex = int;

using SourceArray = std::vector<std::shared_ptr<Source>>;
using ObjectArray = std::vector<std::shared_ptr<Object>>;
using BoundaryArray = std::vector<std::shared_ptr<Boundary>>;
using YeeCellArray = std::vector<std::shared_ptr<YeeCell>>;
using MonitorArray = std::vector<std::shared_ptr<Monitor>>;

enum class PlaneType { XY, YZ, ZX };

inline DoubleArrary3D allocateDoubleArray3D(SpatialIndex nx, SpatialIndex ny,
                                            SpatialIndex nz,
                                            double default_value = 0.0) {
  auto arr{DoubleArrary3D(nx, ny, nz)};
  arr.setConstant(default_value);
  return arr;
}

inline DoubleArrary2D allocateDoubleArray2D(SpatialIndex nx, SpatialIndex ny,
                                            double default_value = 0.0) {
  auto arr{DoubleArrary2D(nx, ny)};
  arr.setConstant(default_value);
  return arr;
}
}  // namespace xfdtd

#endif  // _TYPE_DEFINE_H_
