#include "waveform/gaussian_waveform.h"

#include <cmath>
namespace xfdtd {

GaussianWaveform::GaussianWaveform(double amplitude, double tau, double t_0)
    : _amplitude{amplitude}, _tau{tau}, _t_0{t_0} {};

double GaussianWaveform::getValueByTime(double time) const {
  return exp(-pow((time - _t_0) / _tau, 2)) * _amplitude;
}

}  // namespace xfdtd
