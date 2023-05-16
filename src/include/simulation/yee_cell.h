#ifndef _YEE_CELL_H_
#define _YEE_CELL_H_

#include <memory>

#include "shape/cube.h"

namespace xfdtd {

class YeeCell  {
 public:
  YeeCell(Eigen::Vector3d point, Eigen::Vector3d size, int material_index);
  ~YeeCell() = default;

  inline int getMaterialIndex() const { return _material_index; }

  inline Eigen::Vector3d getCenter() const { return _shape->getCenter(); }

  inline void setMaterialIndex(int material_index) { _material_index = material_index; }

 private:
  std::shared_ptr<Cube> _shape;
  int _material_index{-1};
};

}  // namespace xfdtd

#endif  // _YEE_CELL_H_
#include "shape/cube.h"