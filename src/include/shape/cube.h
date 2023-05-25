#ifndef _CUBE_H_
#define _CUBE_H_

#include <memory>

#include "shape/shape.h"

namespace xfdtd {
/**
 * @brief Axis Aligned Cube
 *
 */
class Cube : public Shape {
 public:
  Cube(PointVector point, PointVector size);
  Cube(const Cube &other) = default;
  Cube(Cube &&other) noexcept = default;
  Cube &operator=(const Cube &other) = default;
  Cube &operator=(Cube &&other) noexcept = default;
  ~Cube() override = default;

  explicit operator std::string() const override;

  std::unique_ptr<Shape> clone() const override;

  bool isPointInside(const PointVector &point) const override;
  std::unique_ptr<Shape> getWrappedBox() const override;

  inline double getXmin() const { return _point(0); }
  inline double getXmax() const { return _point(0) + _size(0); }
  inline double getYmin() const { return _point(1); }
  inline double getYmax() const { return _point(1) + _size(1); }
  inline double getZmin() const { return _point(2); }
  inline double getZmax() const { return _point(2) + _size(2); }

  inline PointVector getCenter() const { return _point + _size / 2; }
  inline PointVector getPoint() const { return _point; }
  inline PointVector getSize() const { return _size; }

 private:
  PointVector _point;
  PointVector _size;
};
}  // namespace xfdtd

#endif  // _CUBE_H_
