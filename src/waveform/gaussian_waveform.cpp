#include "waveform/gaussian_waveform.h"

#include <cmath>
namespace xfdtd {

GaussianWaveform::GaussianWaveform(double amplitude, double tau, double t_0)
    : _amplitude{amplitude}, _tau{tau}, _t_0{t_0} {};

void GaussianWaveform::init(const std::vector<double> &time_array) {
  for (const auto &t : time_array) {
    getAllValues().emplace_back(getValueByTime(t));
  }
}

double GaussianWaveform::getValueByTime(double time) const {
  return exp(-pow((time - _t_0) / _tau, 2)) * _amplitude;
}

}  // namespace xfdtd
