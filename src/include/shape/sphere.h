#ifndef _SPHERE_H_
#define _SPHERE_H_

#include <Eigen/Core>
#include <string_view>

#include "shape/shape.h"
#include "util/type_define.h"

namespace xfdtd {

class Sphere : public Shape {
 public:
  Sphere(PointVector center, double radius);
  Sphere(const Sphere &other) = default;
  Sphere(Sphere &&other) noexcept = default;
  Sphere &operator=(const Sphere &other) = default;
  Sphere &operator=(Sphere &&other) noexcept = default;
  ~Sphere() override = default;

  explicit operator std::string() const override;

  std::unique_ptr<Shape> clone() const override;

  bool isPointInside(const PointVector &point) const override;

  std::unique_ptr<Shape> getWrappedBox() const override;

 private:
  PointVector _center;
  double _radius;
};
}  // namespace xfdtd

#endif  // _SPHERE_H_
