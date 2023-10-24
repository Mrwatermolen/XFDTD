#ifndef _XFDTD_CAPACITOR_H_
#define _XFDTD_CAPACITOR_H_

#include "lumped_element/lumped_element.h"

namespace xfdtd {
class Capacitor : public LumpedElement {
 public:
  Capacitor(std::unique_ptr<Shape> shape, Axis axis, double capacitance);
  ~Capacitor() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<EMF> emf) override;

  void correctFDTDCoff() override;

  void correctE() override;

  void correctH() override;

 private:
  Axis _axis;
  double _capacitance;
  size_t _is, _ie, _js, _je, _ks, _ke;
  double _capacitance_factor;
  xt::xarray<double> _grid_size_a, _grid_size_b, _grid_size_c;
  xt::xarray<double> _da, _db, _dc;
  xt::xarray<double> _beta;

  void correctFDTDCoff(xt::xarray<double>& cece, xt::xarray<double>& cecha,
                       xt::xarray<double>& cechb, const xt::xarray<double>& eps,
                       const xt::xarray<double>& sigma);
};
}  // namespace xfdtd

#endif  // _XFDTD_CAPACITOR_H_
