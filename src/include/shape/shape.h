#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <memory>
#include <string>

#include "util/type_define.h"

namespace xfdtd {

class Shape {
 public:
  virtual ~Shape() = default;

  virtual explicit operator std::string() const;

  virtual std::unique_ptr<Shape> clone() const = 0;

  /**
   * @brief determine if the point is in the shape
   */
  virtual bool isPointInside(const PointVector& point) const = 0;

  /**
   * @brief Get the Wrapped Box object. The wrapped box is a cube.
   *
   * @return std::unique_ptr<Shape>
   */
  virtual std::unique_ptr<Shape> getWrappedBox() const = 0;
};

}  // namespace xfdtd

#endif  // _SHAPE_H_
