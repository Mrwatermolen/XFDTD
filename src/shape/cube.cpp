#include "shape/cube.h"

#include <cstddef>
#include <memory>
#include <utility>

#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"

namespace xfdtd {

Cube::Cube(PointVector point, PointVector size)
    : _point{std::move(point)}, _size{std::move(size)} {};

Cube::operator std::string() const {
  return std::string("Cube: ") + std::to_string(_point(0)) + " " +
         std::to_string(_point(1)) + " " + std::to_string(_point(2)) + " " +
         std::to_string(_size(0)) + " " + std::to_string(_size(1)) + " " +
         std::to_string(_size(2));
}

std::unique_ptr<Shape> Cube::clone() const {
  return std::make_unique<Cube>(*this);
}

bool Cube::isPointInside(const PointVector &point) const {
  return isLessOrEqual(_point(0), point(0), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(_point(1), point(1), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(_point(2), point(2), constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point(0) + _size(0), point(0),
                          constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point(1) + _size(1), point(1),
                          constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(_point(2) + _size(2), point(2),
                          constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Cube::getWrappedBox() const {
  return std::make_unique<Cube>(_point, _size);
}

}  // namespace xfdtd
