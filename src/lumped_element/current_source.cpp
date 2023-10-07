#include "lumped_element/current_source.h"

namespace xfdtd {

CurrentSource::CurrentSource(std::unique_ptr<Cube> cube,
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

CurrentSource::CurrentSource(const CurrentSource& other)
    : LumpedElement{other},
      _orientation{other._orientation},
      _resistance{other._resistance},
      _waveform{other._waveform->clone()} {}

CurrentSource& CurrentSource::operator=(const CurrentSource& other) {
  if (this != &other) {
    LumpedElement::operator=(other);
    _orientation = other._orientation;
    _resistance = other._resistance;
    _waveform = other._waveform->clone();
  }
  return *this;
}

void CurrentSource::init(std::shared_ptr<GridSpace> grid_space,
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
    _resistance_factor = _resistance * (_ke - _ks) * (_je - _js) / (_ie - _is);
    _current_amplitude_factor =
        _waveform->getAmplitude() / ((_ke - _ks) * (_je - _js));
    _da = {fdtd_basic_coff->getDy()};
    _db = {fdtd_basic_coff->getDz()};
    _dc = {fdtd_basic_coff->getDx()};
    if (_orientation == Orientation::XN) {
      _current_amplitude_factor *= -1;
    }
  } else if (_orientation == Orientation::YN ||
             _orientation == Orientation::YP) {
    _resistance_factor = _resistance * (_ke - _ks) * (_ie - _is) / (_je - _js);
    _current_amplitude_factor =
        _waveform->getAmplitude() / ((_ke - _ks) * (_ie - _is));
    _da = {fdtd_basic_coff->getDz()};
    _db = {fdtd_basic_coff->getDx()};
    _dc = {fdtd_basic_coff->getDy()};
    if (_orientation == Orientation::YN) {
      _current_amplitude_factor *= -1;
    }
  } else if (_orientation == Orientation::ZN ||
             _orientation == Orientation::ZP) {
    _resistance_factor = _resistance * (_ie - _is) * (_je - _js) / (_ke - _ks);
    _current_amplitude_factor =
        _waveform->getAmplitude() / ((_ie - _is) * (_je - _js));
    _da = {fdtd_basic_coff->getDx()};
    _db = {fdtd_basic_coff->getDy()};
    _dc = {fdtd_basic_coff->getDz()};
    if (_orientation == Orientation::ZN) {
      _current_amplitude_factor *= -1;
    }
  }

  _beta = (dt * _dc) / (_resistance_factor * _da * _db);
  _waveform->setAmplitude(_current_amplitude_factor);
  _waveform->init(
      xt::linspace(0.5 * dt, (total_time_step - 0.5) * dt, total_time_step));
}

void CurrentSource::correctFDTDCoff() {
  auto fdtd_basic_coff{getFDTDBasicCoff()};

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    correctFDTDCoff(fdtd_basic_coff->getCexe(), fdtd_basic_coff->getCexhy(),
                    fdtd_basic_coff->getCexhz(), fdtd_basic_coff->getEpsX(),
                    fdtd_basic_coff->getSigmaX());

  } else if (_orientation == Orientation::YN ||
             _orientation == Orientation::YP) {
    correctFDTDCoff(fdtd_basic_coff->getCeye(), fdtd_basic_coff->getCeyhz(),
                    fdtd_basic_coff->getCeyhx(), fdtd_basic_coff->getEpsY(),
                    fdtd_basic_coff->getSigmaY());

  } else if (_orientation == Orientation::ZN ||
             _orientation == Orientation::ZP) {
    correctFDTDCoff(fdtd_basic_coff->getCeze(), fdtd_basic_coff->getCezhx(),
                    fdtd_basic_coff->getCezhy(), fdtd_basic_coff->getEpsZ(),
                    fdtd_basic_coff->getSigmaZ());
  }
}

void CurrentSource::correctE() {
  auto current_time_step{getFDTDBasicCoff()->getCurrentTimeStep()};
  auto range_x{xt::range(_is, _ie)};
  auto range_y{xt::range(_js, _je)};
  auto range_z{xt::range(_ks, _ke)};

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    auto ec{xt::view(getEMF()->getEx(), range_x, range_y, range_z)};
    ec += _coff_i * ((*_waveform)(current_time_step));
    return;
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    auto ec{xt::view(getEMF()->getEy(), range_x, range_y, range_z)};
    ec += _coff_i * ((*_waveform)(current_time_step));
    return;
  }

  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    auto ec{xt::view(getEMF()->getEz(), range_x, range_y, range_z)};
    ec += _coff_i * ((*_waveform)(current_time_step));
    return;
  }
}

void CurrentSource::correctH() {}

void CurrentSource::correctFDTDCoff(xt::xarray<double>& cece,
                                    xt::xarray<double>& cecha,
                                    xt::xarray<double>& cechb,
                                    const xt::xarray<double>& eps,
                                    const xt::xarray<double>& sigma) {
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
  _coff_i = -2 * dt / ((2 * eps_view + dt * sigma_view + _beta) * (_db * _da));
}
}  // namespace xfdtd
