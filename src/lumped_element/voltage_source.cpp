#include "lumped_element/voltage_source.h"

#include <memory>
#include <xtensor/xarray.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xview.hpp>

#include "shape/cube.h"
#include "util/type_define.h"
#include "waveform/waveform.h"

namespace xfdtd {

VoltageSource::VoltageSource(std::unique_ptr<Cube> cube,
                             Orientation orientation, double resistance,
                             std::unique_ptr<Waveform> waveform)
    : LumpedElement{std::move(cube)},
      _orientation{orientation},
      _resistance{resistance},
      _waveform{std::move(waveform)} {
  if (_resistance == 0) {
    // avoid  divided by zero
    _resistance = 1e-20;
  }
}

VoltageSource::VoltageSource(const VoltageSource &other)
    : LumpedElement{other},
      _orientation{other._orientation},
      _resistance{other._resistance},
      _waveform{other._waveform->clone()} {}

VoltageSource &VoltageSource::operator=(const VoltageSource &other) {
  if (this != &other) {
    LumpedElement::operator=(other);
    _orientation = other._orientation;
    _resistance = other._resistance;
    _waveform = other._waveform->clone();
  }
  return *this;
}

void VoltageSource::init(std::shared_ptr<GridSpace> grid_space,
                         std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                         std::shared_ptr<EMF> emf) {
  defaultInit(grid_space, fdtd_basic_coff, emf);
  // range: [is, ie), [js, je), [ks, ke)
  auto grid_box{grid_space->getGridBox(getShape())};
  auto dt{fdtd_basic_coff->getDt()};
  auto total_time_step{fdtd_basic_coff->getTotalTimeStep()};

  _is = grid_box->getGridStartIndexX();
  _ie = grid_box->getGridEndIndexX();
  _js = grid_box->getGridStartIndexY();
  _je = grid_box->getGridEndIndexY();
  _ks = grid_box->getGridStartIndexZ();
  _ke = grid_box->getGridEndIndexZ();

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    _resistance_factor = _resistance * (_je - _js) * (_ke - _ks) / (_ie - _is);
    _voltage_amplitude_factor = _waveform->getAmplitude() / (_ie - _is);
    if (_orientation == Orientation::XN) {
      _voltage_amplitude_factor *= -1;
    }
    _da = {fdtd_basic_coff->getDy()};
    _db = {fdtd_basic_coff->getDz()};
    _alpha = _da * _db * _resistance_factor;
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDx() /
            (_resistance_factor * _da * _db);
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YN) {
    _resistance_factor = _resistance * (_ie - _is) * (_ke - _ks) / (_je - _js);
    _voltage_amplitude_factor = _waveform->getAmplitude() / (_je - _js);
    if (_orientation == Orientation::YN) {
      _voltage_amplitude_factor *= -1;
    }
    _da = {fdtd_basic_coff->getDz()};
    _db = {fdtd_basic_coff->getDx()};
    _alpha = _da * _db * _resistance_factor;
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDy() /
            (_resistance_factor * _da * _db);
  }

  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    _resistance_factor = _resistance * (_ie - _is) * (_je - _js) / (_ke - _ks);
    _voltage_amplitude_factor = _waveform->getAmplitude() / (_ke - _ks);
    if (_orientation == Orientation::ZN) {
      _voltage_amplitude_factor *= -1;
    }
    _da = {fdtd_basic_coff->getDx()};
    _db = {fdtd_basic_coff->getDy()};
    _alpha = _da * _db * _resistance_factor;
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDz() / (_alpha);
  }

  _waveform->setAmplitude(_voltage_amplitude_factor);
  _waveform->init(xt::linspace<double>(dt * 0.5, (total_time_step - 0.5) * dt,
                                       total_time_step));
}

void VoltageSource::correctFDTDCoff() {
  auto fdtd_basic_coff{getFDTDBasicCoff()};

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    correctFDTDCoff(fdtd_basic_coff->getCexe(), fdtd_basic_coff->getCexhy(),
                    fdtd_basic_coff->getCexhz(), fdtd_basic_coff->getEpsX(),
                    fdtd_basic_coff->getSigmaX());
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    correctFDTDCoff(fdtd_basic_coff->getCeye(), fdtd_basic_coff->getCeyhz(),
                    fdtd_basic_coff->getCeyhx(), fdtd_basic_coff->getEpsY(),
                    fdtd_basic_coff->getSigmaY());
  }

  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    correctFDTDCoff(fdtd_basic_coff->getCeze(), fdtd_basic_coff->getCezhx(),
                    fdtd_basic_coff->getCezhy(), fdtd_basic_coff->getEpsZ(),
                    fdtd_basic_coff->getSigmaZ());
  }
}

void VoltageSource::correctFDTDCoff(xt::xarray<double> &cece,
                                    xt::xarray<double> &cecha,
                                    xt::xarray<double> &cechb,
                                    const xt::xarray<double> &eps,
                                    const xt::xarray<double> &sigma) {
  auto range_x{xt::range(_is, _ie)};
  auto range_y{xt::range(_js, _je)};
  auto range_z{xt::range(_ks, _ke)};
  auto dt{getFDTDBasicCoff()->getDt()};
  auto cece_view{xt::view(cece, range_x, range_y, range_z)};
  auto cecha_view{xt::view(cecha, range_x, range_y, range_z)};
  auto cechb_view{xt::view(cechb, range_x, range_y, range_z)};
  const auto eps_view{xt::view(eps, range_x, range_y, range_z)};
  const auto sigma_view{xt::view(sigma, range_x, range_y, range_z)};

  cece_view = (2 * eps_view - dt * sigma_view - _beta) /
              (2 * eps_view + dt * sigma_view + _beta);
  cecha_view = -2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * _db);
  cechb_view = 2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * _da);
  _coff_v = -2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * _alpha);
}

void VoltageSource::correctE() {
  auto current_time_step{getFDTDBasicCoff()->getCurrentTimeStep()};
  auto range_x{xt::range(_is, _ie)};
  auto range_y{xt::range(_js, _je)};
  auto range_z{xt::range(_ks, _ke)};

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    auto ec{xt::view(getEMF()->getEx(), range_x, range_y, range_z)};
    ec += _coff_v * ((*_waveform)(current_time_step));
    return;
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    auto ec{xt::view(getEMF()->getEy(), range_x, range_y, range_z)};
    ec += _coff_v * ((*_waveform)(current_time_step));
    return;
  }

  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    auto ec{xt::view(getEMF()->getEz(), range_x, range_y, range_z)};
    ec += _coff_v * ((*_waveform)(current_time_step));
    return;
  }
}

void VoltageSource::correctH() {}
}  // namespace xfdtd
