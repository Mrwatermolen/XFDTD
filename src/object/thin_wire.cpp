#include "object/thin_wire.h"

#include <memory>
#include <xtensor/xview.hpp>

#include "shape/cylinder.h"
#include "util/type_define.h"

namespace xfdtd {

ThinWire::ThinWire(const std::string& name, std::unique_ptr<Cylinder> shape)
    : Object(name, std::move(shape),
             Material::createAir(
                 "placeholder_wire"))  // It doesn't mean that ThinWire is air.
                                       // It just need a material.
{}

void ThinWire::init(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<GridSpace> grid_space,
                    const std::shared_ptr<EMF>& emf) {
  defaultInit(index, fdtd_basic_coff, grid_space, emf);
  if (!grid_space->isUniformGridSpace()) {
    throw std::runtime_error("ThinWire only support uniform grid space");
  }
  auto shape{getShape()};
  auto cylinder{dynamic_cast<Cylinder*>(shape.get())};
  if (cylinder == nullptr) {
    throw std::runtime_error("ThinWire only support Cylinder shape");
  }

  _axis = cylinder->getAxis();
  _radius = cylinder->getRadius();

  if (!isEnoughThin()) {
    throw std::runtime_error("ThinWire is not enough thin");
  }
}

void ThinWire::correctFDTDCoff() {
  // NOTE: Current material type does not support magnetizing.
  auto fdtd_basic_coff{getFDTDBasicCoff()};
  auto grid_space{getGridSpace()};
  auto shape{getShape()};
  auto grid_box{grid_space->getGridBox(shape.get())};
  auto is{grid_box->getGridOriginIndexX()};
  auto ie{grid_box->getGridEndIndexX()};
  auto js{grid_box->getGridOriginIndexY()};
  auto je{grid_box->getGridEndIndexY()};
  auto ks{grid_box->getGridOriginIndexZ()};
  auto ke{grid_box->getGridEndIndexZ()};
  auto x_range{xt::range(is, ie)};
  auto y_range{xt::range(js, je)};
  auto z_range{xt::range(ks, ke)};
  auto dx{grid_space->getGridBaseSizeX()};
  auto dy{grid_space->getGridBaseSizeY()};
  auto dz{grid_space->getGridBaseSizeZ()};
  auto dt = fdtd_basic_coff->getDt();

  if (_axis == Axis::X) {
    // Ex is zero in thin wire
    // js == je
    // ks == ke
    auto cexe_view{xt::view(fdtd_basic_coff->getCexe(), x_range, js, ks)};
    auto cexhy_view{xt::view(fdtd_basic_coff->getCexhy(), x_range, js, ks)};
    auto cexhz_view{xt::view(fdtd_basic_coff->getCexhz(), x_range, js, ks)};
    cexe_view = 0.0;
    cexhy_view = 0.0;
    cexhz_view = 0.0;

    auto khy{dz * std::atan(dy / dy) / dy};
    auto khz{dy * std::atan(dz / dy) / dz};
    auto chyex_view{xt::view(fdtd_basic_coff->getChyex(), x_range, js,
                             xt::range(ks - 1, ks + 1))};
    chyex_view = khy * 2 * chyex_view / std::log(dz / _radius);

    auto chzex_view{xt::view(fdtd_basic_coff->getChzex(), x_range,
                             xt::range(js - 1, js + 1), ks)};
    chzex_view = khz * 2 * chzex_view / std::log(dy / _radius);

    return;
  }

  if (_axis == Axis::Y) {
    // Ey is zero in thin wire
    // is == ie
    // ks == ke
    auto ceye_view{xt::view(fdtd_basic_coff->getCeye(), is, y_range, ks)};
    auto ceyhz_view{xt::view(fdtd_basic_coff->getCeyhz(), is, y_range, ks)};
    auto ceyhx_view{xt::view(fdtd_basic_coff->getCeyhx(), is, y_range, ks)};
    // ceye_view = 0.0;
    // ceyhz_view = 0.0;
    // ceyhx_view = 0.0;

    auto khz{dx * std::atan(dz / dx) / dz};
    auto khx{dz * std::atan(dx / dz) / dx};
    auto chzey_view{xt::view(fdtd_basic_coff->getChzey(),
                             xt::range(is - 1, is + 1), y_range, ks)};
    chzey_view = khz * 2 * chzey_view / std::log(dx / _radius);

    auto chxey_view{xt::view(fdtd_basic_coff->getChxey(), is, y_range,
                             xt::range(ks - 1, ks + 1))};
    chxey_view = khx * 2 * chxey_view / std::log(dz / _radius);

    return;
  }

  if (_axis == Axis::Z) {
    // Ez is zero in thin wire
    // is == ie
    // js == je
    auto ceze_view{xt::view(fdtd_basic_coff->getCeze(), is, js, z_range)};
    auto cezhx_view{xt::view(fdtd_basic_coff->getCezhx(), is, js, z_range)};
    auto cezhy_view{xt::view(fdtd_basic_coff->getCezhy(), is, js, z_range)};
    ceze_view = 0.0;
    cezhx_view = 0.0;
    cezhy_view = 0.0;

    auto khx{dy * std::atan(dx / dy) / dx};
    auto khy{dx * std::atan(dy / dx) / dy};
    auto chxez_view{xt::view(fdtd_basic_coff->getChxez(), is,
                             xt::range(js - 1, js + 1), z_range)};
    chxez_view = khx * 2 * chxez_view / std::log(dy / _radius);

    auto chyez_view{xt::view(fdtd_basic_coff->getChyez(),
                             xt::range(is - 1, is + 1), js, z_range)};
    chyez_view = khy * 2 * chyez_view / std::log(dx / _radius);
    return;
  }
}

bool ThinWire::isEnoughThin() const {
  auto dx{getGridSpace()->getGridBaseSizeX()};
  auto dy{getGridSpace()->getGridBaseSizeY()};
  auto dz{getGridSpace()->getGridBaseSizeZ()};
  if (_axis == Axis::X) {
    return _radius < dy / 2 && _radius < dz / 2;
  }
  if (_axis == Axis::Y) {
    return _radius < dx / 2 && _radius < dz / 2;
  }
  if (_axis == Axis::Z) {
    return _radius < dx / 2 && _radius < dy / 2;
  }
  return false;
}

}  // namespace xfdtd
