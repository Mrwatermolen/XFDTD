
#include "shape/cylinder.h"

#include <cmath>
#include <memory>
#include <utility>

#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {
Cylinder::Cylinder(PointVector center, double raduis, double height)
    : _center{std::move(center)}, _raduis{raduis}, _height{height} {}

std::unique_ptr<Shape> Cylinder::clone() const {
  return std::make_unique<Cylinder>(*this);
}
bool Cylinder::isPointInside(const PointVector& point) const {
  if (!getWrappedBox()->isPointInside(point)) {
    return false;
  }
  auto dis{sqrt(pow(point(0) - _center(0), 2) + pow(point(1) - _center(1), 2))};
  return isLessOrEqual(dis, _raduis, constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Cylinder::getWrappedBox() const {
  return std::make_unique<Cube>(
      PointVector{_center(0) - _raduis, _center(1) - _raduis,
                  _center(2) - _height / 2},
      PointVector{2 * _raduis, 2 * _raduis, _height});
}
}  // namespace xfdtd