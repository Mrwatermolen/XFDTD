#include "shape/cube.h"

#include <cstddef>
#include <memory>
#include <utility>

#include "util/constant.h"
#include "shape/shape.h"
#include "util/float_compare.h"

namespace xfdtd {

Cube::Cube(Eigen::Vector3d point, Eigen::Vector3d size)
    : _point{std::move(point)}, _size{std::move(size)} {};

Cube::operator std::string() const {
  return std::string("Cube: ") + std::to_string(_point.x()) + " " +
         std::to_string(_point.y()) + " " + std::to_string(_point.z()) + " " +
         std::to_string(_size.x()) + " " + std::to_string(_size.y()) + " " +
         std::to_string(_size.z());
}

std::unique_ptr<Shape> Cube::clone() const {
  return std::make_unique<Cube>(*this);
}

bool Cube::isPointInside(const Eigen::Vector3d &point) const {
  return isLessOrEqual(_point.x(), point.x(), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(_point.y(), point.y(), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(_point.z(), point.z(), constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point.x() + _size.x(), point.x(),
                          constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point.x() + _size.y(), point.y(),
                          constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point.x() + _size.z(), point.z(),
                          constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Cube::getWrappedBox() const {
  return std::make_unique<Cube>(_point, _size);
}

}  // namespace xfdtd
