#ifndef _XFDTD_RESISTOR_H_
#define _XFDTD_RESISTOR_H_

#include <memory>

#include "lumped_element/lumped_element.h"
#include "shape/cube.h"
#include "util/type_define.h"

namespace xfdtd {

class Resistor : public LumpedElement {
 public:
  Resistor(std::unique_ptr<Cube> shape, Axis axis, double resistance);
  ~Resistor() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<EMF> emf) override;

  void correctFDTDCoff() override;

  void correctE() override;

  void correctH() override;

 private:
  Axis _axis;
  double _resistance, _resistance_factor;
  xt::xarray<double> _grid_size_a, _grid_size_b, _grid_size_c;
  xt::xarray<double> _da, _db, _dc, _beta;
  size_t _is, _ie, _js, _je, _ks, _ke;

  void correctFDTDCoff(xt::xarray<double> &cece, xt::xarray<double> &cecha,
                       xt::xarray<double> &cechb, const xt::xarray<double> &eps,
                       const xt::xarray<double> &sigma);
};
}  // namespace xfdtd

#endif  // _XFDTD_RESISTOR_H_
