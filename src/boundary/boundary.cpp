#include "boundary/boundary.h"

#include <utility>

namespace xfdtd {

void Boundary::defaultInit(std::shared_ptr<EMF> emf,
                           std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                           std::shared_ptr<GridSpace> grid_space,
                           std::unique_ptr<Shape> shape) {
  _emf = std::move(emf);
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _grid_space = std::move(grid_space);
  _shape = std::move(shape);
}

EMF* Boundary::getEMFInstance() const { return _emf.get(); }

FDTDBasicCoff* Boundary::getFDTDBasicCoff() const {
  return _fdtd_basic_coff.get();
}

GridSpace* Boundary::getGridSpace() const { return _grid_space.get(); }

Shape* Boundary::getShape() const { return _shape.get(); }

}  // namespace xfdtd
