#include "simulation/yee_cell.h"

#include <utility>

namespace xfdtd {
YeeCell::YeeCell(Eigen::Vector3d point, Eigen::Vector3d size,
                 int material_index)
    : _shape{std::make_shared<Cube>(std::move(point), std::move(size))},
      _material_index{material_index} {}
}  // namespace xfdtd