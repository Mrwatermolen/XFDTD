#include "waveform/cosine_modulated_gaussian_waveform.h"

#include <cmath>

#include "util/constant.h"

namespace xfdtd {
CosineModulatedGaussianWaveform::CosineModulatedGaussianWaveform(
    double amplitude, double tau, double t_0, double modulation_frequency)
    : _amplitude{amplitude},
      _tau{tau},
      _t_0{t_0},
      _modulation_frequency{modulation_frequency} {}

double CosineModulatedGaussianWaveform::getValueByTime(double time) const {
  return _amplitude *
         std::cos(2 * constant::PI * _modulation_frequency * (time - _t_0)) *
         std::exp(-pow((time - _t_0) / _tau, 2));
}
}  // namespace xfdtd
