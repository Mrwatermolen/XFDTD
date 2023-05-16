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
  Cube(Eigen::Vector3d point, Eigen::Vector3d size);
  Cube(const Cube &other) = default;
  Cube(Cube &&other) noexcept = default;
  Cube &operator=(const Cube &other) = default;
  Cube &operator=(Cube &&other) noexcept = default;
  ~Cube() override = default;

  explicit operator std::string() const override;

  std::unique_ptr<Shape> clone() const override;

  bool isPointInside(const Eigen::Vector3d &point) const override;
  std::unique_ptr<Shape> getWrappedBox() const override;

  inline double getXmin() const { return _point.x(); }
  inline double getXmax() const { return _point.x() + _size.x(); }
  inline double getYmin() const { return _point.y(); }
  inline double getYmax() const { return _point.y() + _size.y(); }
  inline double getZmin() const { return _point.z(); }
  inline double getZmax() const { return _point.z() + _size.z(); }

  inline Eigen::Vector3d getCenter() const { return _point + _size / 2; }
  inline Eigen::Vector3d getPoint() const { return _point; }
  inline Eigen::Vector3d getSize() const { return _size; }

 private:
  Eigen::Vector3d _point;
  Eigen::Vector3d _size;
};
}  // namespace xfdtd

#endif  // _CUBE_H_
