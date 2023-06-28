#include "shape/sphere.h"

#include <cmath>
#include <string>
#include <utility>
#include <xtensor/xnorm.hpp>

#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"

namespace xfdtd {
Sphere::Sphere(PointVector center, double radius)
    : _center{std::move(center)}, _radius{radius} {}

Sphere::operator std::string() const {
  return std::string("Sphere: ") + std::to_string(_center(0)) + " " +
         std::to_string(_center(1)) + " " + std::to_string(_center(2)) + " " +
         std::to_string(_radius);
}

std::unique_ptr<Shape> Sphere::clone() const {
  return std::make_unique<Sphere>(*this);
}

bool Sphere::isPointInside(const PointVector &point) const {
  auto dist{xt::eval(xt::norm_l2(point - _center))()};
  return isLessOrEqual(dist, _radius, constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Sphere::getWrappedBox() const {
  return std::make_unique<Cube>(
      (_center - PointVector{_radius / 2, _radius / 2, _radius / 2}),
      PointVector{_radius, _radius, _radius});
}
}  // namespace xfdtd
