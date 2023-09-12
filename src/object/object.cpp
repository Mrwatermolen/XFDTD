#include "object/object.h"

#include <memory>
#include <xtensor/xindex_view.hpp>

#include "material/material.h"
#include "shape/shape.h"
#include "util/type_define.h"

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

void Object::correctCece(xt::xarray<bool>& mask, EFTA& cece, double dt) {
  auto eps = _material->getPermittivityE();
  auto sigma_e = _material->getElectricalConductivity();
  xt::filter(cece, mask) = (2 * eps - dt * sigma_e) / (2 * eps + dt * sigma_e);
}

void Object::correctCecha(xt::xarray<bool>& mask, EFTA& cecha, double db,
                          double dt) {
  auto eps = _material->getPermittivityE();
  auto sigma_e = _material->getElectricalConductivity();
  xt::filter(cecha, mask) = -(2 * dt / db) / (2 * eps + dt * sigma_e);
}

void Object::correctCechb(xt::xarray<bool>& mask, EFTA& ceahb, double da,
                          double dt) {
  auto eps = _material->getPermittivityE();
  auto sigma_e = _material->getElectricalConductivity();
  xt::filter(ceahb, mask) = (2 * dt / da) / (2 * eps + dt * sigma_e);
}

void Object::correctChch(xt::xarray<bool>& mask, EFTA& chch, double dt) {
  auto mu{_material->getPermeabilityM()};
  auto sigma_m{_material->getMagneticConductivity()};
  xt::filter(chch, mask) = (2 * mu - dt * sigma_m) / (2 * mu + dt * sigma_m);
}

void Object::correctChcea(xt::xarray<bool>& mask, EFTA& chaha, double db,
                          double dt) {
  auto mu{_material->getPermeabilityM()};
  auto sigma_m{_material->getMagneticConductivity()};
  xt::filter(chaha, mask) = (2 * dt / db) / (2 * mu + dt * sigma_m);
}

void Object::correctChceb(xt::xarray<bool>& mask, EFTA& chahb, double da,
                          double dt) {
  auto mu{_material->getPermeabilityM()};
  auto sigma_m{_material->getMagneticConductivity()};
  xt::filter(chahb, mask) = -(2 * dt / da) / (2 * mu + dt * sigma_m);
}

}  // namespace xfdtd
