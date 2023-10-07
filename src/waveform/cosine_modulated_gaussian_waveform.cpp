#include "waveform/cosine_modulated_gaussian_waveform.h"

#include <cmath>

#include "util/constant.h"
#include "waveform/waveform.h"

namespace xfdtd {
CosineModulatedGaussianWaveform::CosineModulatedGaussianWaveform(
    double amplitude, double tau, double t_0, double modulation_frequency)
    : Waveform{amplitude},
      _tau{tau},
      _t_0{t_0},
      _modulation_frequency{modulation_frequency} {}

CosineModulatedGaussianWaveform::CosineModulatedGaussianWaveform(
    double tau, double t_0, double modulation_frequency)
    : _tau{tau}, _t_0{t_0}, _modulation_frequency{modulation_frequency} {}

std::unique_ptr<Waveform> CosineModulatedGaussianWaveform::clone() const {
  return std::make_unique<CosineModulatedGaussianWaveform>(*this);
}

double CosineModulatedGaussianWaveform::getValueByTime(double time) const {
  return getAmplitude() *
         std::cos(2 * constant::PI * _modulation_frequency * (time - _t_0)) *
         std::exp(-pow((time - _t_0) / _tau, 2));
}

void CosineModulatedGaussianWaveform::init(xt::xarray<double> time_array) {
  setValue(
      getAmplitude() *
      xt::cos(2 * constant::PI * _modulation_frequency * (time_array - _t_0)) *
      xt::exp(-xt::pow((time_array - _t_0) / _tau, 2)));
}
}  // namespace xfdtd
