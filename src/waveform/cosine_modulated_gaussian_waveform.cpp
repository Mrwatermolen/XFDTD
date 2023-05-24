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

void CosineModulatedGaussianWaveform::init(
    const std::vector<double> &time_array) {
  for (const auto &t : time_array) {
    getAllValues().emplace_back(getValueByTime(t));
  }
}

double CosineModulatedGaussianWaveform::getValueByTime(double time) const {
  return _amplitude *
         std::cos(2 * constant::PI * _modulation_frequency * time) *
         std::exp(-4 * constant::PI * pow((time - _t_0) / _tau, 2));
}
}  // namespace xfdtd
