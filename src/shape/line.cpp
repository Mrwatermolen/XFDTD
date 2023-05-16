#include "shape/line.h"

#include <Eigen/Dense>
#include <iostream>
#include <memory>

#include "util/constant.h"
#include "shape/cube.h"
#include "util/float_compare.h"

namespace xfdtd {
Line::Line(Eigen::Vector3d start, Eigen::Vector3d end)
    : _start(std::move(start)), _end(std::move(end)) {}

Line::operator std::string() const {
  return "Line: " + std::to_string(_start.x()) + " " +
         std::to_string(_start.y()) + " " + std::to_string(_start.z()) +
         " -> " + std::to_string(_end.x()) + " " + std::to_string(_end.y()) +
         " " + std::to_string(_end.z());
}

std::unique_ptr<Shape> Line::clone() const {
  return std::make_unique<Line>(*this);
}

double Line::getNorm() const { return (_end - _start).norm(); }

bool Line::isPointInside(const Eigen::Vector3d &point) const {
  return isGreaterOrEqual(point.x(), _start.x(), constant::TOLERABLE_EPSILON) &&
         isGreaterOrEqual(point.y(), _start.y(), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(point.x(), _end.x(), constant::TOLERABLE_EPSILON) &&
         isLessOrEqual(point.y(), _end.y(), constant::TOLERABLE_EPSILON);
}

bool Line::isVertical(const Line &line) const {
  return isEqual((_end - _start).dot((line._end - line._start)), 0,
                 constant::TOLERABLE_EPSILON);
}

bool Line::isParallel(const Line &line) const {
  return isEqual((_end - _start).cross((line._end - line._start)).norm(), 0,
                 constant::TOLERABLE_EPSILON);
}

std::unique_ptr<Shape> Line::getWrappedBox() const {
  return std::make_unique<Cube>(_start, _end);
}
}  // namespace xfdtd
