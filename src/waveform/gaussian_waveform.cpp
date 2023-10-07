#include "waveform/gaussian_waveform.h"

#include <cmath>

#include "waveform/waveform.h"
namespace xfdtd {

GaussianWaveform::GaussianWaveform(double amplitude, double tau, double t_0)
    : Waveform{amplitude}, _tau{tau}, _t_0{t_0} {};

GaussianWaveform::GaussianWaveform(double tau, double t_0)
    : _tau{tau}, _t_0{t_0} {};

std::unique_ptr<Waveform> GaussianWaveform::clone() const {
  return std::make_unique<GaussianWaveform>(*this);
}

double GaussianWaveform::getValueByTime(double time) const {
  return getAmplitude() * std::exp(-std::pow((time - _t_0) / _tau, 2));
}

void GaussianWaveform::init(xt::xarray<double> time_array) {
  setValue(getAmplitude() * xt::exp(-xt::pow((time_array - _t_0) / _tau, 2)));
}

}  // namespace xfdtd
