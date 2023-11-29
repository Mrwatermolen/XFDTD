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
           double phi_inc, double psi, std::shared_ptr<Waveform> waveform)
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

void TFSF::defaultInitTFSF(std::shared_ptr<const GridSpace> grid_space,
                           std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                           std::shared_ptr<EMF> emf) {
  if (!grid_space->isUniformGridSpace()) {
    throw std::runtime_error("TFSF doesn't support non-uniform grid spacing.");
  }

  _grid_space = std::move(grid_space);
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _emf = std::move(emf);
  _tfsf_grid_box =
      std::make_unique<GridBox>(_distance_x, _distance_y, _distance_z,
                                _grid_space->getGridNumX() - 2 * _distance_x,
                                _grid_space->getGridNumY() - 2 * _distance_y,
                                _grid_space->getGridNumZ() - 2 * _distance_z);

  auto time_arr{xt::arange<double>(0, getTotalTimeStep()) * getDt()};
  initIncidentWaveForm(time_arr);
}

xt::xarray<double> TFSF::getIncidentWaveValue() const {
  return _e_0 * _waveform->getValue();
}

void TFSF::initIncidentWaveForm(xt::xarray<double> time) {
  _waveform->init(std::move(time));
}

xt::xarray<std::complex<double>> TFSF::getIncidentWaveFourierTransform(
    const xt::xarray<double> &frequencies) const {
  return dft(_waveform->getValue(), getDt(), frequencies);
}

std::tuple<xt::xarray<double>, xt::xarray<std::complex<double>>>
TFSF::getIncidentWaveFastFourierTransform() const {
  auto temp{_waveform->getValue()};
  auto n{2 * getTotalTimeStep() - 1};
  auto additive_shape{n - temp.size()};
  auto w{xt::concatenate(
      xt::xtuple(std::move(temp), xt::zeros<double>({additive_shape})), 0)};
  auto v{xt::fftw::fftshift(xt::fftw::fft(w))};
  auto frequencies = xt::fftw::rfftfreq(w.size(), getDt());
  return {frequencies, v};
}

double TFSF::getIncidentWaveValueByTimeStep(size_t t) const {
  return _e_0 * (*_waveform)(t);
}

}  // namespace xfdtd