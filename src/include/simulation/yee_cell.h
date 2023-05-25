#ifndef _YEE_CELL_H_
#define _YEE_CELL_H_

#include <memory>
#include <tuple>

#include "shape/cube.h"
#include "util/type_define.h"

namespace xfdtd {

class YeeCell {
 public:
  YeeCell(PointVector point, PointVector size, int material_index,
          SpatialIndex x = -1, SpatialIndex y = -1, SpatialIndex z = -1);
  ~YeeCell() = default;

  inline int getMaterialIndex() const { return _material_index; }

  inline PointVector getCenter() const { return _shape->getCenter(); }

  inline void setMaterialIndex(int material_index) {
    _material_index = material_index;
  }

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getGridXYZIndex()
      const {
    return std::make_tuple(_x, _y, _z);
  }

 private:
  std::shared_ptr<Cube> _shape;
  int _material_index{-1};
  SpatialIndex _x;
  SpatialIndex _y;
  SpatialIndex _z;
};

}  // namespace xfdtd

#endif  // _YEE_CELL_H_
#include "shape/cube.h"
