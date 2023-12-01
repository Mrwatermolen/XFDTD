#ifndef _XFDTD_FLUX_H_
#define _XFDTD_FLUX_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_box.h"
#include "grid/grid_space.h"
#include "shape/shape.h"

namespace xfdtd {

class Flux {
 public:
  Flux(std::unique_ptr<Shape> shape, EMComponent component);

  Flux(std::unique_ptr<GridBox> grid_box, EMComponent component);

  virtual ~Flux() = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<const EMF> emf);

  void update();

  void output();

  xt::xarray<double>& getValue() { return _value; }

 private:
  std::unique_ptr<Shape> _shape;
  EMComponent _component;
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<const EMF> _emf;
  std::unique_ptr<GridBox> _grid_box;
  xt::xarray<double> _value;

  void captureX();
  void captureY();
  void captureZ();
};

}  // namespace xfdtd

#endif  // _XFDTD_FLUX_H_
