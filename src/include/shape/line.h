#ifndef __LINE_H__
#define __LINE_H__

#include "shape/curve.h"

namespace xfdtd {
class Line : public Curve {
 public:
  Line(PointVector start, PointVector end);
  Line(const Line &line) = default;
  Line(Line &&line) noexcept = default;
  Line &operator=(const Line &line) = default;
  Line &operator=(Line &&line) = default;
  ~Line() override = default;

  explicit operator std::string() const override;

  std::unique_ptr<Shape> clone() const override;

  double getNorm() const;

  bool isPointInside(const PointVector &point) const override;
  bool isVertical(const Line &line) const;
  bool isParallel(const Line &line) const;

  std::unique_ptr<Shape> getWrappedBox() const override;

 private:
  PointVector _start;
  PointVector _end;
};
}  //   namespace xfdtd

#endif  // __LINE_H__
