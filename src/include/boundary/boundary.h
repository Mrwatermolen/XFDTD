#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include <memory>
#include <utility>

#include "electromagnetic_field/electromagnetic_field.h"
#include "util/type_define.h"

namespace xfdtd {
class Simulation;
enum class Orientation { XN, XP, YN, YP, ZN, ZP };
class Boundary {
 public:
  Boundary() = default;
  Boundary(const Boundary &) = default;
  Boundary(Boundary &&) noexcept = default;
  Boundary &operator=(const Boundary &) = default;
  Boundary &operator=(Boundary &&) noexcept = default;
  virtual ~Boundary() = default;

  virtual SpatialIndex getSize() const = 0;
  virtual Orientation getOrientation() const = 0;

  inline std::shared_ptr<EMF> getEMFInstance() const {
    if (_emf == nullptr) {
      throw std::runtime_error{"EMF is not set"};
    }
    return _emf;
  }

  void setEMFInstance(std::shared_ptr<EMF> emf) { _emf = std::move(emf); }

  virtual void init(Simulation *simulation) = 0;
  virtual void updateH() = 0;
  virtual void updateE() = 0;

 private:
  std::shared_ptr<EMF> _emf;
};
}  // namespace xfdtd

#endif  // _BOUNDARY_H_
