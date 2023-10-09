
#include "shape/cylinder.h"

#include <cmath>
#include <limits>
#include <memory>
#include <utility>

#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {
Cylinder::Cylinder(Axis axis, PointVector center, double radius, double height)
    : _axis{axis},
      _center{std::move(center)},
      _radius{radius},
      _height{height} {}

std::unique_ptr<Shape> Cylinder::clone() const {
  return std::make_unique<Cylinder>(*this);
}

bool Cylinder::isPointInside(const PointVector& point) const {
  if (!getWrappedBox()->isPointInside(point)) {
    return false;
  }
  double dis{std::numeric_limits<double>::max()};
  switch (_axis) {
    case Axis::X:
      dis = std::sqrt(std::pow(point(1) - _center(1), 2) +
                      std::pow(point(2) - _center(2), 2));
      break;
    case Axis::Y:
      dis = std::sqrt(std::pow(point(0) - _center(0), 2) +
                      std::pow(point(2) - _center(2), 2));
      break;
    case Axis::Z:
      dis = std::sqrt(std::pow(point(0) - _center(0), 2) +
                      std::pow(point(1) - _center(1), 2));
      break;
  }
  return isLessOrEqual(dis, _radius, constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Cylinder::getWrappedBox() const {
  switch (_axis) {
    case Axis::X:
      return std::make_unique<Cube>(
          PointVector{_center(0) - _height / 2, _center(1) - _radius,
                      _center(2) - _radius},
          PointVector{_height, 2 * _radius, 2 * _radius});
    case Axis::Y:
      return std::make_unique<Cube>(
          PointVector{_center(0) - _radius, _center(1) - _height / 2,
                      _center(2) - _radius},
          PointVector{2 * _radius, _height, 2 * _radius});
    case Axis::Z:
      return std::make_unique<Cube>(
          PointVector{_center(0) - _radius, _center(1) - _radius,
                      _center(2) - _height / 2},
          PointVector{2 * _radius, 2 * _radius, _height});
  }
}

Axis Cylinder::getAxis() const { return _axis; }

double Cylinder::getRadius() const { return _radius; }

double Cylinder::getHeight() const { return _height; }

PointVector Cylinder::getCenter() const { return _center; }

}  // namespace xfdtd
