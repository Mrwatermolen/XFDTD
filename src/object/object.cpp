#include "object/object.h"

#include <memory>

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
  _material = other._material->clone();
}

Object& Object::operator=(const Object& other) {
  if (this != &other) {
    _name = other._name;
    _shape = other._shape->clone();
    _material = other._material->clone();
  }
  return *this;
}

Object::operator std::string() const {
  return std::string("Object: ") + _name + "\n" +
         static_cast<std::string>(*_shape) + "\n" +
         static_cast<std::string>(*_material);
}

std::unique_ptr<Object> Object::clone() const {
  return std::make_unique<Object>(_name, _shape->clone(), _material->clone());
}

}  // namespace xfdtd
