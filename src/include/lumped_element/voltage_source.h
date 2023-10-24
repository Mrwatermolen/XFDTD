#ifndef _XFDTD_VOLTAGE_SOURCE_H_
#define _XFDTD_VOLTAGE_SOURCE_H_

#include "lumped_element/lumped_element.h"
#include "waveform/waveform.h"

namespace xfdtd {

class VoltageSource : public LumpedElement {
 public:
  /**
   * @brief Construct a new Voltage Source:: Voltage Source object
   *
   * @param cube shape of voltage source
   * @param orientation
   * @param resistance the internal resistance of voltage source. If it is zero,
   * it means this is a hard voltage source.
   * @param waveform the waveform of voltage source
   */

  VoltageSource(std::unique_ptr<Cube> cube, Orientation orientation,
                double resistance, std::unique_ptr<Waveform> waveform);
  VoltageSource(const VoltageSource& other);
  VoltageSource& operator=(const VoltageSource& other);
  VoltageSource(VoltageSource&& other) noexcept = default;
  VoltageSource& operator=(VoltageSource&& other) noexcept = default;
  ~VoltageSource() override = default;

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
  size_t _is, _ie, _js, _je, _ks, _ke;
  double _resistance_factor;
  double _voltage_amplitude_factor;
  xt::xarray<double> _grid_size_a, _grid_size_b, _grid_size_c;
  xt::xarray<double> _da, _db, _dc;
  xt::xarray<double> _alpha;
  xt::xarray<double> _beta;
  xt::xarray<double> _coff_v;

  void correctFDTDCoff(xt::xarray<double>& cece, xt::xarray<double>& cecha,
                       xt::xarray<double>& cechb, const xt::xarray<double>& eps,
                       const xt::xarray<double>& sigma);
};

}  // namespace xfdtd

#endif  // _XFDTD_VOLTAGE_SOURCE_H_