#include "lumped_element/inductor.h"

#include "util/type_define.h"

namespace xfdtd {

Inductor::Inductor(std::unique_ptr<Cube> shape, Axis axis, double inductance)
    : LumpedElement{std::move(shape)}, _axis{axis}, _inductance{inductance} {}

void Inductor::init(std::shared_ptr<GridSpace> grid_space,
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
    _inductance_factor = _inductance * (_je - _js) * (_ke - _ks) / (_ie - _is);
    _da = {fdtd_basic_coff->getDy()};
    _db = {fdtd_basic_coff->getDz()};
    _dc = {fdtd_basic_coff->getDx()};
  }
  if (_axis == Axis::Y) {
    _inductance_factor = _inductance * (_ie - _is) * (_ke - _ks) / (_je - _js);
    _da = {fdtd_basic_coff->getDz()};
    _db = {fdtd_basic_coff->getDx()};
    _dc = {fdtd_basic_coff->getDy()};
  }
  if (_axis == Axis::Z) {
    _inductance_factor = _inductance * (_ie - _is) * (_je - _js) / (_ke - _ks);
    _da = {fdtd_basic_coff->getDx()};
    _db = {fdtd_basic_coff->getDy()};
    _dc = {fdtd_basic_coff->getDz()};
  }
  auto dt{getFDTDBasicCoff()->getDt()};
  _cjcec = dt * _dc / (_inductance_factor * _da * _db);
}

void Inductor::correctFDTDCoff() {
  auto x_range = xt::range(_is, _ie);
  auto y_range = xt::range(_js, _je);
  auto z_range = xt::range(_ks, _ke);
  auto dt{getFDTDBasicCoff()->getDt()};
  auto nx{getGridSpace()->getGridNumX()};
  auto ny{getGridSpace()->getGridNumY()};
  auto nz{getGridSpace()->getGridNumZ()};
  if (_axis == Axis::X) {
    const auto eps_view{
        xt::view(getFDTDBasicCoff()->getEpsX(), x_range, y_range, z_range)};
    const auto sigma_view{
        xt::view(getFDTDBasicCoff()->getSigmaX(), x_range, y_range, z_range)};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    getEMF()->allocateJx(nx, ny + 1, nz + 1);
    return;
  }

  if (_axis == Axis::Y) {
    const auto eps_view{
        xt::view(getFDTDBasicCoff()->getEpsY(), x_range, y_range, z_range)};
    const auto sigma_view{
        xt::view(getFDTDBasicCoff()->getSigmaY(), x_range, y_range, z_range)};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    getEMF()->allocateJy(nx + 1, ny, nz + 1);
    return;
  }

  if (_axis == Axis::Z) {
    const auto eps_view{
        xt::view(getFDTDBasicCoff()->getEpsZ(), x_range, y_range, z_range)};
    const auto sigma_view{
        xt::view(getFDTDBasicCoff()->getSigmaZ(), x_range, y_range, z_range)};
    _cecjc = -(2 * dt) / (2 * eps_view + dt * sigma_view);
    getEMF()->allocateJz(nx + 1, ny + 1, nz);
    return;
  }
}

void Inductor::correctE() {
  if (_axis == Axis::X) {
    auto& jx{getEMF()->getJx()};
    auto& ex{getEMF()->getEx()};
    auto jx_view{xt::view(jx, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    auto ex_view{xt::view(ex, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    ex_view += _cecjc * jx_view;
    jx_view = (_cjcec * ex_view + jx_view);
    return;
  }

  if (_axis == Axis::Y) {
    auto& jy{getEMF()->getJy()};
    auto& ey{getEMF()->getEy()};
    auto jy_view{xt::view(jy, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    auto ey_view{xt::view(ey, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    ey_view += _cecjc * jy_view;
    jy_view = (_cjcec * ey_view + jy_view);
    return;
  }

  if (_axis == Axis::Z) {
    auto& jz{getEMF()->getJz()};
    auto& ez{getEMF()->getEz()};
    auto jz_view{xt::view(jz, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    auto ez_view{xt::view(ez, xt::range(_is, _ie), xt::range(_js, _je),
                          xt::range(_ks, _ke))};
    ez_view += _cecjc * jz_view;
    jz_view = (_cjcec * ez_view + jz_view);
    return;
  }
}

void Inductor::correctH() {}

}  // namespace xfdtd
