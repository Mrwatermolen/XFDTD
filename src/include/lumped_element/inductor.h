#ifndef _XFDTD_INDUCTOR_H_
#define _XFDTD_INDUCTOR_H_

#include "lumped_element/lumped_element.h"
#include "shape/cube.h"

namespace xfdtd {

class Inductor : public LumpedElement {
 public:
  Inductor(std::unique_ptr<Cube> shape, Axis axis, double inductance);
  ~Inductor() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<EMF> emf) override;

  void correctFDTDCoff() override;

  void correctE() override;

  void correctH() override;

 private:
  Axis _axis;
  double _da, _db, _dc;
  size_t _is, _ie, _js, _je, _ks, _ke;
  double _inductance;
  double _inductance_factor;
  double _cjcec;
  xt::xarray<double> _cecjc;
};

}  // namespace xfdtd

#endif  // _XFDTD_INDUCTOR_H_
