#include "waveform/gaussian_waveform.h"

namespace xfdtd {

GaussianWaveform::GaussianWaveform(double amplitude, double tau, double t_0)
    : _amplitude{amplitude}, _tau{tau}, _t_0{t_0} {};

void GaussianWaveform::init(const std::vector<double> &time_array) {
  _size = time_array.size();
  _values = std::make_unique<double[]>(_size);
  for (size_t i = 0; i < _size; ++i) {
    auto temp = (time_array[i] - _t_0) / _tau;
    auto temp1 = pow(temp, 2);
    auto temp2 = exp(-temp1);
    _values[i] = getValueByTime(time_array[i]);
  }
}

double GaussianWaveform::getValueByTime(double time) const {
  return exp(-pow((time - _t_0) / _tau, 2)) * _amplitude;
}

std::unique_ptr<Waveform> GaussianWaveform::clone() const {
  return std::make_unique<GaussianWaveform>(*this);
}

}  // namespace xfdtd
