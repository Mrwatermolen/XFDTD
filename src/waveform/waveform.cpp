#include "waveform/waveform.h"

namespace xfdtd {

const std::vector<double>& Waveform::getAllValues() const { return _values; }

std::vector<double>& Waveform::getAllValues() { return _values; }

double Waveform::getValue(size_t time_step) const { return _values[time_step]; }
}  // namespace xfdtd
