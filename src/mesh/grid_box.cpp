#include "mesh/grid_box.h"

namespace xfdtd {
GridBox::GridBox(SpatialIndex start_x, SpatialIndex start_y,
                 SpatialIndex start_z, SpatialIndex nx, SpatialIndex ny,
                 SpatialIndex nz)
    : _start_x{start_x},
      _start_y{start_y},
      _start_z{start_z},
      _nx{nx},
      _ny{ny},
      _nz{nz} {}
}  // namespace xfdtd