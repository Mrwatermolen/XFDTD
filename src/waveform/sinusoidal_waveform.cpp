#include "waveform/sinusoidal_waveform.h"

#include <cmath>

#include "util/constant.h"

namespace xfdtd {

SinusoidalWaveform::SinusoidalWaveform(double amplitude, double frequency,
                                       double phase)
    : _amplitude(amplitude), _frequency(frequency), _phase(phase) {}

double SinusoidalWaveform::getValueByTime(double time) const {
  return _amplitude * std::sin(2 * constant::PI * _frequency * time + _phase);
}

}  // namespace xfdtd
