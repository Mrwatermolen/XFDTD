#include "shape/sphere.h"

#include <cmath>
#include <string>
#include <utility>

#include "util/constant.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/float_compare.h"

namespace xfdtd {
Sphere::Sphere(Eigen::Vector3d center, double radius)
    : _center{std::move(center)}, _radius{radius} {}

Sphere::operator std::string() const {
  return std::string("Sphere: ") + std::to_string(_center.x()) + " " +
         std::to_string(_center.y()) + " " + std::to_string(_center.z()) + " " +
         std::to_string(_radius);
}

std::unique_ptr<Shape> Sphere::clone() const {
  return std::make_unique<Sphere>(*this);
}

bool Sphere::isPointInside(const Eigen::Vector3d &point) const {
  auto dist{(point - _center).norm()};
  return isLessOrEqual(dist, _radius, constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Sphere::getWrappedBox() const {
  return std::make_unique<Cube>(
      (_center - Eigen::Vector3d{_radius / 2, _radius / 2, _radius / 2}),
      Eigen::Vector3d{_radius, _radius, _radius});
}
}  // namespace xfdtd
