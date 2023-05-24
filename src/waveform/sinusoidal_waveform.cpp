#include "waveform/sinusoidal_waveform.h"

#include <vector>

#include "util/constant.h"

namespace xfdtd {

SinusoidalWaveform::SinusoidalWaveform(double amplitude, double frequency,
                                       double phase)
    : _amplitude(amplitude), _frequency(frequency), _phase(phase) {}

void SinusoidalWaveform::init(const std::vector<double>& time_array) {
  for (const auto& t : time_array) {
    getAllValues().emplace_back(getValueByTime(t));
  }
}

double SinusoidalWaveform::getValueByTime(double time) const {
  return _amplitude * std::sin(2 * constant::PI * _frequency * time + _phase);
}

}  // namespace xfdtd
