#include "waveform/sinusoidal_waveform.h"

#include <vector>

#include "util/constant.h"

namespace xfdtd {

SinusoidalWaveform::SinusoidalWaveform(double amplitude, double frequency,
                                       double phase)
    : _amplitude(amplitude), _frequency(frequency), _phase(phase) {}

void SinusoidalWaveform::init(const std::vector<double>& time_array) {
  _size = time_array.size();
  _values = std::make_unique<double[]>(_size);
  for (size_t i{0}; i < _size; ++i) {
    _values[i] = getValueByTime(time_array[i]);
  }
}

double SinusoidalWaveform::getValueByTime(double time) const {
  return _amplitude * std::sin(2 * constant::PI * _frequency * time + _phase);
}

std::unique_ptr<Waveform> SinusoidalWaveform::clone() const {
  return std::make_unique<SinusoidalWaveform>(*this);
}

}  // namespace xfdtd
