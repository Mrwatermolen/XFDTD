#ifndef _GRID_BOX_H_
#define _GRID_BOX_H_

#include "util/type_define.h"

namespace xfdtd {
class GridBox {
 public:
  GridBox(SpatialIndex start_x, SpatialIndex start_y, SpatialIndex start_z,
          SpatialIndex nx, SpatialIndex ny, SpatialIndex nz);
  ~GridBox() = default;

  inline SpatialIndex getStartIndexX() const { return _start_x; }
  inline SpatialIndex getEndIndexX() const { return _start_x + _nx; }
  inline SpatialIndex getStartIndexY() const { return _start_y; }
  inline SpatialIndex getEndIndexY() const { return _start_y + _ny; }
  inline SpatialIndex getStartIndexZ() const { return _start_z; }
  inline SpatialIndex getEndIndexZ() const { return _start_z + _nz; }
  inline SpatialIndex getNx() const { return _nx; }
  inline SpatialIndex getNy() const { return _ny; }
  inline SpatialIndex getNz() const { return _nz; }
  inline SpatialIndex getCenterIndexX() const { return _start_x + _nx / 2; }
  inline SpatialIndex getCenterIndexY() const { return _start_y + _ny / 2; }
  inline SpatialIndex getCenterIndexZ() const { return _start_z + _nz / 2; }

 private:
  SpatialIndex _start_x;
  SpatialIndex _start_y;
  SpatialIndex _start_z;
  SpatialIndex _nx;
  SpatialIndex _ny;
  SpatialIndex _nz;
};
}  // namespace xfdtd

#endif  // _GRID_BOX_H_