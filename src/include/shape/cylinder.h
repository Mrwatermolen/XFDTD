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
  Cylinder(Axis axis, PointVector center, double radius, double height);
  ~Cylinder() override = default;

  std::unique_ptr<Shape> clone() const override;

  bool isPointInside(const PointVector& point) const override;

  std::unique_ptr<Shape> getWrappedBox() const override;

  Axis getAxis() const;

  double getRadius() const;

  double getHeight() const;

  PointVector getCenter() const;

 private:
  Axis _axis;
  PointVector _center;
  double _radius;
  double _height;
};
}  // namespace xfdtd

#endif  // _CYLINDER_H_