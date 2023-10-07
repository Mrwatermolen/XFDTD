#include "lumped_element/resistor.h"

#include "util/type_define.h"

namespace xfdtd {

Resistor::Resistor(std::unique_ptr<Cube> shape, Axis axis, double resistance)
    : LumpedElement{std::move(shape)}, _axis{axis}, _resistance{resistance} {}

void Resistor::init(std::shared_ptr<GridSpace> grid_space,
                    std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<EMF> emf) {
  defaultInit(grid_space, fdtd_basic_coff, emf);
  auto grid_box{grid_space->getGridBox(getShape())};
  auto dt{fdtd_basic_coff->getDt()};
  auto total_time_step{fdtd_basic_coff->getTotalTimeStep()};

  _is = grid_box->getGridStartIndexX();
  _ie = grid_box->getGridEndIndexX();
  _js = grid_box->getGridStartIndexY();
  _je = grid_box->getGridEndIndexY();
  _ks = grid_box->getGridStartIndexZ();
  _ke = grid_box->getGridEndIndexZ();

  if (_axis == Axis::X) {
    _resistance_factor = _resistance * (_je - _js) * (_ke - _ks) / (_ie - _is);
    _da = {fdtd_basic_coff->getDy()};
    _db = {fdtd_basic_coff->getDz()};
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDx() /
            (_resistance_factor * _da * _db);
  }

  if (_axis == Axis::Y) {
    _resistance_factor = _resistance * (_ie - _is) * (_ke - _ks) / (_je - _js);
    _da = {fdtd_basic_coff->getDz()};
    _db = {fdtd_basic_coff->getDx()};
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDy() /
            (_resistance_factor * _da * _db);
  }

  if (_axis == Axis::Z) {
    _resistance_factor = _resistance * (_ie - _is) * (_je - _js) / (_ke - _ks);
    _da = {fdtd_basic_coff->getDx()};
    _db = {fdtd_basic_coff->getDy()};
    _beta = fdtd_basic_coff->getDt() * fdtd_basic_coff->getDz() /
            (_resistance_factor * _da * _db);
  }
}

void Resistor::correctFDTDCoff() {
  auto fdtd_basic_coff{getFDTDBasicCoff()};
  if (_axis == Axis::X) {
    correctFDTDCoff(fdtd_basic_coff->getCexe(), fdtd_basic_coff->getCexhy(),
                    fdtd_basic_coff->getCexhz(), fdtd_basic_coff->getEpsX(),
                    fdtd_basic_coff->getSigmaX());
  }

  if (_axis == Axis::Y) {
    correctFDTDCoff(fdtd_basic_coff->getCeye(), fdtd_basic_coff->getCeyhz(),
                    fdtd_basic_coff->getCeyhx(), fdtd_basic_coff->getEpsY(),
                    fdtd_basic_coff->getSigmaY());
  }

  if (_axis == Axis::Z) {
    correctFDTDCoff(fdtd_basic_coff->getCeze(), fdtd_basic_coff->getCezhx(),
                    fdtd_basic_coff->getCezhy(), fdtd_basic_coff->getEpsZ(),
                    fdtd_basic_coff->getSigmaZ());
  }
}

void Resistor::correctFDTDCoff(xt::xarray<double> &cece,
                               xt::xarray<double> &cecha,
                               xt::xarray<double> &cechb,
                               const xt::xarray<double> &eps,
                               const xt::xarray<double> &sigma) {
  auto range_x = xt::range(_is, _ie);
  auto range_y = xt::range(_js, _je);
  auto range_z = xt::range(_ks, _ke);
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
}

void Resistor::correctE(){};

void Resistor::correctH(){};

}  // namespace xfdtd
