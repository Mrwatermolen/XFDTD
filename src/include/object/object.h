#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <memory>
#include <string_view>
#include <tuple>

#include "material/dispersive_material.h"
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
  virtual ~Object() = default;

  explicit operator std::string() const;

  std::unique_ptr<Object> clone() const;

  inline bool isPointInside(const PointVector& point) {
    return _shape->isPointInside(point);
  }

  inline std::unique_ptr<Shape> getWrappedBox() const {
    return _shape->getWrappedBox();
  }

  inline bool isDispersion() { return _material->isDispersion(); }

  /**
   * @brief get the electromagnetic properties of the Material
   * @return a tuple of (epsilon, mu, sigma, conductivity)
   */
  inline std::tuple<double, double, double, double>
  getElectromagneticProperties() {
    return _material->getElectromagneticProperties();
  }

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf) {
    if (isDispersion()) {
      auto m{_material.get()};
      m->init(dt, dl, emf);
    }
  }

  virtual void correctCece(xt::xarray<bool>& mask, EFTA& cece, double dt);
  virtual void correctCecha(xt::xarray<bool>& mask, EFTA& cecha, double db,
                            double dt);
  virtual void correctCechb(xt::xarray<bool>& mask, EFTA& ceahb, double da,
                            double dt);
  virtual void correctChch(xt::xarray<bool>& mask, EFTA& chch, double dt);
  virtual void correctChcea(xt::xarray<bool>& mask, EFTA& chaha, double db,
                            double dt);
  virtual void correctChceb(xt::xarray<bool>& mask, EFTA& chahb, double da,
                            double dt);

  void updateEx(int i, int j, int k) { _material->updateEx(i, j, k); }
  void updateEy(int i, int j, int k) { _material->updateEy(i, j, k); }
  void updateEz(int i, int j, int k) { _material->updateEz(i, j, k); }

 private:
  std::string _name;
  std::unique_ptr<Shape> _shape;
  std::unique_ptr<Material> _material;
};

}  // namespace xfdtd

#endif  // _OBJECT_H_
