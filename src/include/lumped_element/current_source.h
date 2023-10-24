#ifndef _XFDTD_CURRENT_SOURCE_H_
#define _XFDTD_CURRENT_SOURCE_H_

#include "lumped_element/lumped_element.h"
#include "waveform/waveform.h"

namespace xfdtd {

class CurrentSource : public LumpedElement {
 public:
  CurrentSource(std::unique_ptr<Cube> cube, Orientation orientation,
                double resistance, std::unique_ptr<Waveform> waveform);
  CurrentSource(const CurrentSource& other);
  CurrentSource& operator=(const CurrentSource& other);
  CurrentSource(CurrentSource&& other) noexcept = default;
  CurrentSource& operator=(CurrentSource&& other) noexcept = default;
  ~CurrentSource() override = default;

  void init(std::shared_ptr<GridSpace> grid_space,
            std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<EMF> emf) override;

  void correctFDTDCoff() override;

  void correctE() override;

  void correctH() override;

 private:
  Orientation _orientation;
  double _resistance;
  std::unique_ptr<Waveform> _waveform;
  // double _da, _db, _dc;
  size_t _is, _ie, _js, _je, _ks, _ke;
  double _resistance_factor;
  double _current_amplitude_factor;
  // double _beta;
  xt::xarray<double> _grid_size_a, _grid_size_b, _grid_size_c;
  xt::xarray<double> _da, _db, _dc;
  xt::xarray<double> _beta;
  xt::xarray<double> _coff_i;

  void correctFDTDCoff(xt::xarray<double>& cece, xt::xarray<double>& cecha,
                       xt::xarray<double>& cechb, const xt::xarray<double>& eps,
                       const xt::xarray<double>& sigma);
};

}  // namespace xfdtd

#endif  // _XFDTD_CURRENT_SOURCE_H_
