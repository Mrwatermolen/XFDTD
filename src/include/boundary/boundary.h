#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

namespace xfdtd {
enum class Orientation { XN, XP, YN, YP, ZN, ZP };
class Boundary {
 public:
  Boundary() = default;
  Boundary(const Boundary &) = default;
  Boundary(Boundary &&) noexcept = default;
  Boundary &operator=(const Boundary &) = default;
  Boundary &operator=(Boundary &&) noexcept = default;
  virtual ~Boundary() = default;

  virtual int getSize() const = 0;
  virtual Orientation getOrientation() const = 0;
  // virtual void updateE() = 0;
  // virtual void updateH() = 0;
};
}  // namespace xfdtd

#endif  // _BOUNDARY_H_
