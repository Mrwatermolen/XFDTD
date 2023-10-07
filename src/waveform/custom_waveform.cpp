#include "waveform/custom_waveform.h"

namespace xfdtd {

CustomWaveform::CustomWaveform(double amplitude,
                               std::function<double(double)> f)
    : Waveform(amplitude), _f(std::move(f)) {}

std::unique_ptr<Waveform> CustomWaveform::clone() const {
  return std::make_unique<CustomWaveform>(*this);
}

void CustomWaveform::init(xt::xarray<double> time_array) {
  xt::xarray<double> val{xt::zeros_like(time_array)};
  for (auto i = 0; i < time_array.size(); ++i) {
    val(i) = getAmplitude() * _f(time_array(i));
  }
  setValue(val);
}

double CustomWaveform::getValueByTime(double time) const {
  return getAmplitude() * _f(time);
}

}  // namespace xfdtd
