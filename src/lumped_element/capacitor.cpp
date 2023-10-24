#include "lumped_element/capacitor.h"

#include "util/type_define.h"

namespace xfdtd {

Capacitor::Capacitor(std::unique_ptr<Shape> shape, Axis axis,
                     double capacitance)
    : LumpedElement{std::move(shape)}, _axis{axis}, _capacitance{capacitance} {}

void Capacitor::init(std::shared_ptr<GridSpace> grid_space,
                     std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                     std::shared_ptr<EMF> emf) {
  defaultInit(grid_space, fdtd_basic_coff, emf);
  // range: [is, ie), [js, je), [ks, ke)
  auto grid_box{grid_space->getGridBox(getShape())};

  _is = grid_box->getGridStartIndexX();
  _ie = grid_box->getGridEndIndexX();
  _js = grid_box->getGridStartIndexY();
  _je = grid_box->getGridEndIndexY();
  _ks = grid_box->getGridStartIndexZ();
  _ke = grid_box->getGridEndIndexZ();

  if (_axis == Axis::X) {
    _capacitance_factor =
        _capacitance * (_ie - _is) / ((_je - _js) * (_ke - _ks));
    _grid_size_a =
        xt::view(grid_space->getGridSizeArrayHY(), xt::range(_js, _je));
    _grid_size_b =
        xt::view(grid_space->getGridSizeArrayHZ(), xt::range(_ks, _ke));
    _grid_size_c =
        xt::view(grid_space->getGridSizeArrayEX(), xt::range(_is, _ie));

    auto [dx, dy, dz] = xt::meshgrid(_grid_size_c, _grid_size_a, _grid_size_b);
    _da = dy;
    _db = dz;
    _dc = dx;
  }

  if (_axis == Axis::Y) {
    _capacitance_factor =
        _capacitance * (_je - _js) / ((_ie - _is) * (_ke - _ks));
    _grid_size_a =
        xt::view(grid_space->getGridSizeArrayHZ(), xt::range(_ks, _ke));
    _grid_size_b =
        xt::view(grid_space->getGridSizeArrayHX(), xt::range(_is, _ie));
    _grid_size_c =
        xt::view(grid_space->getGridSizeArrayEY(), xt::range(_js, _je));

    auto [dx, dy, dz] = xt::meshgrid(_grid_size_b, _grid_size_c, _grid_size_a);
    _da = dz;
    _da = dx;
    _dc = dy;
  }

  if (_axis == Axis::Z) {
    _capacitance_factor =
        _capacitance * (_ke - _ks) / ((_ie - _is) * (_je - _js));
    _grid_size_a =
        xt::view(grid_space->getGridSizeArrayHX(), xt::range(_is, _ie));
    _grid_size_b =
        xt::view(grid_space->getGridSizeArrayHY(), xt::range(_js, _je));
    _grid_size_c =
        xt::view(grid_space->getGridSizeArrayEZ(), xt::range(_ks, _ke));

    auto [dx, dy, dz] = xt::meshgrid(_grid_size_a, _grid_size_b, _grid_size_c);
    _da = dx;
    _db = dy;
    _dc = dz;
  }
  _beta = 2 * _capacitance_factor * _dc / (_da * _db);
}

void Capacitor::correctFDTDCoff() {
  if (_axis == Axis::X) {
    correctFDTDCoff(
        getFDTDBasicCoff()->getCexe(), getFDTDBasicCoff()->getCexhy(),
        getFDTDBasicCoff()->getCexhz(), getFDTDBasicCoff()->getEpsX(),
        getFDTDBasicCoff()->getSigmaX());
  }
  if (_axis == Axis::Y) {
    correctFDTDCoff(
        getFDTDBasicCoff()->getCeye(), getFDTDBasicCoff()->getCeyhz(),
        getFDTDBasicCoff()->getCeyhx(), getFDTDBasicCoff()->getEpsY(),
        getFDTDBasicCoff()->getSigmaY());
  }
  if (_axis == Axis::Z) {
    correctFDTDCoff(
        getFDTDBasicCoff()->getCeze(), getFDTDBasicCoff()->getCezhx(),
        getFDTDBasicCoff()->getCezhy(), getFDTDBasicCoff()->getEpsZ(),
        getFDTDBasicCoff()->getSigmaZ());
  }
}

void Capacitor::correctE() {}

void Capacitor::correctH() {}

void Capacitor::correctFDTDCoff(xt::xarray<double>& cece,
                                xt::xarray<double>& cecha,
                                xt::xarray<double>& cechb,
                                const xt::xarray<double>& eps,
                                const xt::xarray<double>& sigma) {
  auto x_range = xt::range(_is, _ie);
  auto y_range = xt::range(_js, _je);
  auto z_range = xt::range(_ks, _ke);
  auto dt{getFDTDBasicCoff()->getDt()};

  auto cece_view = xt::view(cece, x_range, y_range, z_range);
  auto cecha_view = xt::view(cecha, x_range, y_range, z_range);
  auto cechb_view = xt::view(cechb, x_range, y_range, z_range);
  const auto eps_view = xt::view(eps, x_range, y_range, z_range);
  const auto sigma_view = xt::view(sigma, x_range, y_range, z_range);

  cece_view = (2 * eps_view - dt * sigma_view + _beta) /
              (2 * eps_view + dt * sigma_view + _beta);
  cecha_view = -(2 * dt) / ((2 * eps_view + dt * sigma_view + _beta) * _db);
  cechb_view = (2 * dt) / ((2 * eps_view + dt * sigma_view + _beta) * _da);
}

}  // namespace xfdtd
