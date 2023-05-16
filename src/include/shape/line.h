#ifndef __LINE_H__
#define __LINE_H__

#include "shape/curve.h"

namespace xfdtd {
class Line : public Curve {
 public:
  Line(Eigen::Vector3d start, Eigen::Vector3d end);
  Line(const Line &line) = default;
  Line(Line &&line) noexcept = default;
  Line &operator=(const Line &line) = default;
  Line &operator=(Line &&line) = default;
  ~Line() override = default;

  explicit operator std::string() const override;

  std::unique_ptr<Shape> clone() const override;

  double getNorm() const;

  bool isPointInside(const Eigen::Vector3d &point) const override;
  bool isVertical(const Line &line) const;
  bool isParallel(const Line &line) const;

  std::unique_ptr<Shape> getWrappedBox() const override;

 private:
  Eigen::Vector3d _start;
  Eigen::Vector3d _end;
};
}  //   namespace xfdtd

#endif  // __LINE_H__
