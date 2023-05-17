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
}  // namespace xfdtd

#endif  // _TYPE_DEFINE_H_
