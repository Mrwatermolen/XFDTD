#include "tfsf/tfsf.h"

#include <cmath>
#include <complex>
#include <memory>
#include <tuple>
#include <utility>

#include "util/dft.h"
#include "util/type_define.h"
#include "xtensor-fftw/basic_double.hpp"
#include "xtensor-fftw/helper.hpp"

namespace xfdtd {
TFSF::TFSF(SpatialIndex distance_x, SpatialIndex distance_y,
           SpatialIndex distance_z, double e_0, double theta_inc,
           double phi_inc, double psi, std::unique_ptr<Waveform> waveform)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _e_0{e_0},
      _theta_inc{theta_inc},
      _phi_inc{phi_inc},
      _psi{psi},
      _sin_theta_inc{std::sin(theta_inc)},
      _cos_theta_inc{std::cos(theta_inc)},
      _sin_phi_inc{std::sin(phi_inc)},
      _cos_phi_inc{std::cos(phi_inc)},
      _sin_psi{std::sin(psi)},
      _cos_psi{std::cos(psi)},
      _k{PointVector{_sin_theta_inc * _cos_phi_inc,
                     _sin_theta_inc * _sin_phi_inc, _cos_theta_inc}},
      _waveform{std::move(waveform)} {}

void TFSF::defaultInitTFSF(const GridSpace *grid_space,
                           const FDTDBasicCoff *fdtd_basic_coff,
                           std::shared_ptr<EMF> emf,
                           std::unique_ptr<GridBox> tfsf_grid_box) {
  _dx = grid_space->getGridBaseSizeX();
  _dy = grid_space->getGridBaseSizeY();
  _dz = grid_space->getGridBaseSizeZ();
  _dt = fdtd_basic_coff->getDt();
  _total_time_steps = fdtd_basic_coff->getTotalTimeStep();
  if (_dx != _dy || _dx != _dz) {
    throw std::runtime_error("TFSF doesn't support non-uniform grid spacing.");
  }
  _emf = std::move(emf);
  _tfsf_grid_box = std::move(tfsf_grid_box);
}

xt::xarray<double> TFSF::getIncidentWaveValue() const {
  return _e_0 * _waveform->getValue();
}

void TFSF::initIncidentWaveForm(xt::xarray<double> time) {
  _waveform->init(std::move(time));
}

xt::xarray<std::complex<double>> TFSF::getIncidentWaveFourierTransform(
    const xt::xarray<double> &frequencies) const {
  return dft(_waveform->getValue(), _dt, frequencies);
}

std::tuple<xt::xarray<double>, xt::xarray<std::complex<double>>>
TFSF::getIncidentWaveFastFourierTransform() const {
  auto temp{_waveform->getValue()};
  auto n{2 * _total_time_steps - 1};
  auto additive_shape{n - temp.size()};
  auto w{xt::concatenate(
      xt::xtuple(std::move(temp), xt::zeros<double>({additive_shape})), 0)};
  auto v{xt::fftw::fftshift(xt::fftw::fft(w))};
  auto frequencies = xt::fftw::rfftfreq(w.size(), _dt);
  return {frequencies, v};
}

double TFSF::getIncidentWaveValueByTimeStep(size_t t) const {
  return _e_0 * (*_waveform)(t);
}

}  // namespace xfdtd