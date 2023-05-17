#include "simulation/yee_cell.h"

#include <utility>

namespace xfdtd {
YeeCell::YeeCell(Eigen::Vector3d point, Eigen::Vector3d size,
                 int material_index, SpatialIndex x, SpatialIndex y,
                 SpatialIndex z)
    : _shape{std::make_shared<Cube>(std::move(point), std::move(size))},
      _material_index{material_index},
      _x{x},
      _y{y},
      _z{z} {}
}  // namespace xfdtd
