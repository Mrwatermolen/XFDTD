#include "waveform/sinusoidal_waveform.h"

#include <cmath>

#include "util/constant.h"
#include "waveform/waveform.h"

namespace xfdtd {

SinusoidalWaveform::SinusoidalWaveform(double amplitude, double frequency,
                                       double phase)
    : Waveform(amplitude), _frequency(frequency), _phase(phase) {}

SinusoidalWaveform::SinusoidalWaveform(double frequency, double phase)
    : _frequency(frequency), _phase(phase) {}

std::unique_ptr<Waveform> SinusoidalWaveform::clone() const {
  return std::make_unique<SinusoidalWaveform>(*this);
}

void SinusoidalWaveform::init(xt::xarray<double> time_array) {
  setValue(getAmplitude() *
           xt::sin(2 * constant::PI * _frequency * time_array + _phase));
}

double SinusoidalWaveform::getValueByTime(double time) const {
  return getAmplitude() *
         std::sin(2 * constant::PI * _frequency * time + _phase);
}

}  // namespace xfdtd
