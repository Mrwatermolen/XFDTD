#ifndef _TYPE_DEFINE_H_
#define _TYPE_DEFINE_H_

#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "object/object.h"
#include "simulation/yee_cell.h"

namespace xfdtd {

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
using MaterialArray = std::vector<std::shared_ptr<Object>>;
using BoundaryArray = std::vector<std::shared_ptr<Boundary>>;
using YeeCellArray = std::vector<std::shared_ptr<YeeCell>>;

}  // namespace xfdtd

#endif  // _TYPE_DEFINE_H_