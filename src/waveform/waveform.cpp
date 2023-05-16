#include "waveform/waveform.h"

namespace xfdtd {

std::shared_ptr<double[]> Waveform::getAllValues() const { return _values; }

std::shared_ptr<double[]> Waveform::getValues(size_t start, size_t end) const {
  std::shared_ptr<double[]> values{new double[end - start]};
  for (size_t i = start; i < end; ++i) {
    values[i - start] = _values[i];
  }
  return values;
}

double Waveform::getValue(size_t time_step) { return _values[time_step]; }
}  // namespace xfdtd
