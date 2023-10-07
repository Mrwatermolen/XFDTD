#ifndef _XFDTD_LUMPED_ELEMENT_H_
#define _XFDTD_LUMPED_ELEMENT_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
#include "shape/shape.h"

namespace xfdtd {

class LumpedElement {
 public:
  explicit LumpedElement(std::unique_ptr<Shape> shape);
  LumpedElement(const LumpedElement& other);
  LumpedElement& operator=(const LumpedElement&);
  LumpedElement(LumpedElement&&) noexcept = default;
  LumpedElement& operator=(LumpedElement&&) noexcept = default;
  virtual ~LumpedElement() = default;

  virtual void init(std::shared_ptr<GridSpace> grid_space,
                    std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<EMF> emf) = 0;

  virtual void correctFDTDCoff() = 0;

  virtual void correctE() = 0;

  virtual void correctH() = 0;

 protected:
  void defaultInit(std::shared_ptr<GridSpace> grid_space,
                   std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                   std::shared_ptr<EMF> emf);

  const Shape* getShape() const;

  GridSpace* getGridSpace() const;

  FDTDBasicCoff* getFDTDBasicCoff() const;

  EMF* getEMF() const;

 private:
  std::unique_ptr<Shape> _shape;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<EMF> _emf;
};

}  // namespace xfdtd

#endif  // _XFDTD_LUMPED_ELEMENT_H_
