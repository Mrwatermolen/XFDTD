#ifndef _YEE_CELL_H_
#define _YEE_CELL_H_

#include <tuple>

#include "util/type_define.h"

namespace xfdtd {

class YeeCell {
 public:
  YeeCell() = default;
  YeeCell(SpatialIndex x, SpatialIndex y, SpatialIndex z, double dl,
          double min_x, double min_y, double min_z, int material_index = -1);
  ~YeeCell() = default;

  inline int getMaterialIndex() const { return _material_index; }

  PointVector getCenter() const;

  inline void setMaterialIndex(int material_index) {
    _material_index = material_index;
  }

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getGridXYZIndex()
      const {
    return std::make_tuple(_x, _y, _z);
  }

 private:
  int _material_index{-1};
  SpatialIndex _x;
  SpatialIndex _y;
  SpatialIndex _z;
  double _dl;
  double _min_x, _min_y, _min_z;
};

}  // namespace xfdtd

#endif  // _YEE_CELL_H_
