#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <memory>
#include <string_view>
#include <tuple>

#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
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

  virtual bool isPointInside(const PointVector& point) {
    return _shape->isPointInside(point);
  }

  inline std::unique_ptr<Shape> getWrappedBox() const {
    return _shape->getWrappedBox();
  }

  std::unique_ptr<Shape> getShape() const { return _shape->clone(); }

  inline bool isDispersion() { return _material->isDispersion(); }

  /**
   * @brief get the electromagnetic properties of the Material
   * @return a tuple of (epsilon, mu, sigma, conductivity)
   */
  inline std::tuple<double, double, double, double>
  getElectromagneticProperties() {
    return _material->getElectromagneticProperties();
  }

  virtual void init(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<GridSpace> grid_space,
                    const std::shared_ptr<EMF>& emf);

  virtual void correctFDTDCoff();

  void updateEx(size_t i, size_t j, size_t k);

  void updateEy(size_t i, size_t j, size_t k);

  void updateEz(size_t i, size_t j, size_t k);

  void updateHx(size_t i, size_t j, size_t k);

  void updateHy(size_t i, size_t j, size_t k);

  void updateHz(size_t i, size_t j, size_t k);

 protected:
  void defaultInit(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                   std::shared_ptr<GridSpace> grid_space,
                   const std::shared_ptr<EMF>& emf);

  std::shared_ptr<FDTDBasicCoff> getFDTDBasicCoff() const;

  std::shared_ptr<GridSpace> getGridSpace() const;

  Shape* getShapeRawPoint() const;

 private:
  std::string _name;
  std::unique_ptr<Shape> _shape;
  std::unique_ptr<Material> _material;
  std::shared_ptr<FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<GridSpace> _grid_space;
  double _dt;
  double _dl;
};

}  // namespace xfdtd

#endif  // _OBJECT_H_
