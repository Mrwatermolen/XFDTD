#include "lumped_element/lumped_element.h"

#include <memory>
#include <utility>

namespace xfdtd {

LumpedElement::LumpedElement(std::unique_ptr<Shape> shape)
    : _shape{std::move(shape)} {}

LumpedElement::LumpedElement(const LumpedElement& other) {
  if (&other == this) {
    return;
  }
  _shape = other._shape->clone();
  _grid_space = other._grid_space;
  _fdtd_basic_coff = other._fdtd_basic_coff;
  _emf = other._emf;
}

LumpedElement& LumpedElement::operator=(const LumpedElement& other) {
  if (&other == this) {
    return *this;
  }
  _shape = other._shape->clone();
  _grid_space = other._grid_space;
  _fdtd_basic_coff = other._fdtd_basic_coff;
  _emf = other._emf;
  return *this;
}

void LumpedElement::defaultInit(std::shared_ptr<GridSpace> grid_space,
                                std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                                std::shared_ptr<EMF> emf) {
  _grid_space = std::move(grid_space);
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _emf = std::move(emf);
}

const Shape* LumpedElement::getShape() const { return _shape.get(); }

GridSpace* LumpedElement::getGridSpace() const { return _grid_space.get(); }

FDTDBasicCoff* LumpedElement::getFDTDBasicCoff() const {
  return _fdtd_basic_coff.get();
}

EMF* LumpedElement::getEMF() const { return _emf.get(); }
}  // namespace xfdtd
