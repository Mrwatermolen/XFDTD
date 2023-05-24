#ifndef _CYLINDER_H_
#define _CYLINDER_H_

#include <memory>

#include "shape/shape.h"
#include "util/type_define.h"
namespace xfdtd {

/**
 * @brief only support that the cylinder is parallel to the z axis
 *
 */
class Cylinder : public Shape {
 public:
  Cylinder(PointVector center, double raduis, double height);
  ~Cylinder() override = default;

  std::unique_ptr<Shape> clone() const override;

  bool isPointInside(const Eigen::Vector3d& point) const override;

  std::unique_ptr<Shape> getWrappedBox() const override;

 private:
  PointVector _center;
  double _raduis;
  double _height;
};
}  // namespace xfdtd

#endif  // _CYLINDER_H_