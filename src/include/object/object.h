#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <memory>
#include <string_view>
#include <tuple>

#include "material/material.h"
#include "shape/shape.h"

namespace xfdtd {
class Object {
 public:
  Object(std::string_view name, std::unique_ptr<Shape> shape,
         std::unique_ptr<Material> material);
  Object(std::string_view name, std::unique_ptr<Shape> shape,
         Material material);
  Object(const Object& other);
  Object(Object&& other) noexcept = default;
  Object& operator=(const Object& other);
  Object& operator=(Object&& other) noexcept = default;

  explicit operator std::string() const;

  std::unique_ptr<Object> clone() const;

  inline bool isPointInside(const PointVector& ponit) {
    return _shape->isPointInside(ponit);
  }

  inline std::unique_ptr<Shape> getWrappedBox() const {
    return _shape->getWrappedBox();
  }

  /**
   * @brief get the electromagnetic properties of the Material
   * @return a tuple of (epsilon, mu, sigma, conductivity)
   */
  inline std::tuple<double, double, double, double>
  getElectromagneticProperties() {
    return _material->getElectromagneticProperties();
  }

 private:
  std::string _name;
  std::unique_ptr<Shape> _shape;
  std::unique_ptr<Material> _material;
};

}  // namespace xfdtd

#endif  // _OBJECT_H_
