#include "object/object.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xnpy.hpp>

#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "material/material.h"
#include "shape/shape.h"

namespace xfdtd {

Object::Object(std::string_view name, std::unique_ptr<Shape> shape,
               std::unique_ptr<Material> material)
    : _name{name}, _shape{std::move(shape)}, _material{std::move(material)} {}

Object::Object(std::string_view name, std::unique_ptr<Shape> shape,
               Material material)
    : _name{name},
      _shape{std::move(shape)},
      _material{std::make_unique<Material>(std::move(material))} {}

Object::Object(const Object& other) {
  _name = other._name;
  _shape = other._shape->clone();
  _material = std::make_unique<Material>(*other._material);
}

Object& Object::operator=(const Object& other) {
  if (this != &other) {
    _name = other._name;
    _shape = other._shape->clone();
    _material = std::make_unique<Material>(*other._material);
  }
  return *this;
}

Object::operator std::string() const {
  return std::string("Object: ") + _name + "\n" +
         static_cast<std::string>(*_shape) + "\n" +
         static_cast<std::string>(*_material);
}

std::unique_ptr<Object> Object::clone() const {
  return std::make_unique<Object>(_name, _shape->clone(), *_material);
}

void Object::init(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                  std::shared_ptr<GridSpace> grid_space,
                  const std::shared_ptr<EMF>& emf) {
  defaultInit(index, std::move(fdtd_basic_coff), std::move(grid_space), emf);

  auto grid_box = _grid_space->getGridBox(_shape.get());
  auto [eps, mu, sigma, sigma_m] = _material->getElectromagneticProperties();
  auto mask = _grid_space->getShapeMask(_shape.get());
  xt::filter(_fdtd_basic_coff->getEpsX(), mask) = eps;
  xt::filter(_fdtd_basic_coff->getEpsY(), mask) = eps;
  xt::filter(_fdtd_basic_coff->getEpsZ(), mask) = eps;
  xt::filter(_fdtd_basic_coff->getMuX(), mask) = mu;
  xt::filter(_fdtd_basic_coff->getMuY(), mask) = mu;
  xt::filter(_fdtd_basic_coff->getMuZ(), mask) = mu;
  xt::filter(_fdtd_basic_coff->getSigmaX(), mask) = sigma;
  xt::filter(_fdtd_basic_coff->getSigmaY(), mask) = sigma;
  xt::filter(_fdtd_basic_coff->getSigmaZ(), mask) = sigma;
  xt::filter(_fdtd_basic_coff->getSigmaMx(), mask) = sigma_m;
  xt::filter(_fdtd_basic_coff->getSigmaMy(), mask) = sigma_m;
  xt::filter(_fdtd_basic_coff->getSigmaMz(), mask) = sigma_m;
}

void Object::correctFDTDCoff() {}

void Object::updateEx(size_t i, size_t j, size_t k) {
  _material->updateEx(i, j, k);
}

void Object::updateEy(size_t i, size_t j, size_t k) {
  _material->updateEy(i, j, k);
}

void Object::updateEz(size_t i, size_t j, size_t k) {
  _material->updateEz(i, j, k);
}

void Object::updateHx(size_t i, size_t j, size_t k) {
  _material->updateHx(i, j, k);
}

void Object::updateHy(size_t i, size_t j, size_t k) {
  _material->updateHy(i, j, k);
}

void Object::updateHz(size_t i, size_t j, size_t k) {
  _material->updateHz(i, j, k);
}

void Object::defaultInit(int index,
                         std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                         std::shared_ptr<GridSpace> grid_space,
                         const std::shared_ptr<EMF>& emf) {
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _grid_space = std::move(grid_space);
  _dt = _fdtd_basic_coff->getDt();

  if (isDispersion()) {
    auto m{_material.get()};
    m->init(_grid_space, _fdtd_basic_coff, emf);
    auto grid{_grid_space->getGridView(_shape.get())};
    std::for_each(grid.begin(), grid.end(),
                  [index](auto&& e) { e->setMaterialIndex(index); });
  }
}

std::shared_ptr<FDTDBasicCoff> Object::getFDTDBasicCoff() const {
  return _fdtd_basic_coff;
}

std::shared_ptr<GridSpace> Object::getGridSpace() const { return _grid_space; }

Shape* Object::getShapeRawPoint() const { return _shape.get(); }

}  // namespace xfdtd
