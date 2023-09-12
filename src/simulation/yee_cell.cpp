#include "simulation/yee_cell.h"

namespace xfdtd {
YeeCell::YeeCell(SpatialIndex x, SpatialIndex y, SpatialIndex z, double dl,
                 double min_x, double min_y, double min_z, int material_index)
    : _material_index{material_index},
      _x{x},
      _y{y},
      _z{z},
      _dl{dl},
      _min_x{min_x},
      _min_y{min_y},
      _min_z{min_z} {}

PointVector YeeCell::getCenter() const {
  if (_z == -1) {
    return PointVector{_min_x + _x * _dl + _dl / 2, _min_y + _y * _dl + _dl / 2,
                       0};
  }
  return PointVector{_min_x + _x * _dl + _dl / 2, _min_y + _y * _dl + _dl / 2,
                     _min_z + _z * _dl + _dl / 2};
}
}  // namespace xfdtd
