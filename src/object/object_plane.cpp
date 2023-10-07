#include "object/object_plane.h"

namespace xfdtd {
ObjectPlane::ObjectPlane(std::string_view name, std::unique_ptr<Shape> shape,
                         std::unique_ptr<Material> material)
    : Object{name, std::move(shape), std::move(material)} {}

void ObjectPlane::init(int index,
                       std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                       std::shared_ptr<GridSpace> grid_space,
                       const std::shared_ptr<EMF>& emf) {
  defaultInit(index, fdtd_basic_coff, grid_space, emf);
  auto grid_box = getGridSpace()->getGridBox(getShapeRawPoint());
  auto [eps, mu, sigma, sigma_m] = getElectromagneticProperties();
  if (grid_box->getGridOriginIndexX() == grid_box->getGridEndIndexX()) {
    // plane
    auto is = grid_box->getGridOriginIndexX();
    auto ie = grid_box->getGridEndIndexX();
    auto js = grid_box->getGridOriginIndexY();
    auto je = grid_box->getGridEndIndexY();
    auto ks = grid_box->getGridOriginIndexZ();
    auto ke = grid_box->getGridEndIndexZ();
    xt::view(getFDTDBasicCoff()->getSigmaY(), is, xt::range(js, je),
             xt::range(ks, ke + 1)) = sigma;
    xt::view(getFDTDBasicCoff()->getSigmaZ(), is, xt::range(js, je + 1),
             xt::range(ks, ke)) = sigma;
    return;
  }
  if (grid_box->getGridOriginIndexY() == grid_box->getGridEndIndexY()) {
    auto is = grid_box->getGridOriginIndexX();
    auto ie = grid_box->getGridEndIndexX();
    auto js = grid_box->getGridOriginIndexY();
    auto je = grid_box->getGridEndIndexY();
    auto ks = grid_box->getGridOriginIndexZ();
    auto ke = grid_box->getGridEndIndexZ();
    xt::view(getFDTDBasicCoff()->getSigmaX(), xt::range(is, ie), js,
             xt::range(ks, ke + 1)) = sigma;
    xt::view(getFDTDBasicCoff()->getSigmaZ(), xt::range(is, ie + 1), js,
             xt::range(ks, ke)) = sigma;
    return;
  }
  if (grid_box->getGridOriginIndexZ() == grid_box->getGridEndIndexZ()) {
    auto is = grid_box->getGridOriginIndexX();
    auto ie = grid_box->getGridEndIndexX();
    auto js = grid_box->getGridOriginIndexY();
    auto je = grid_box->getGridEndIndexY();
    auto ks = grid_box->getGridOriginIndexZ();
    auto ke = grid_box->getGridEndIndexZ();
    xt::view(getFDTDBasicCoff()->getSigmaX(), xt::range(is, ie),
             xt::range(js, je + 1), ks) = sigma;
    xt::view(getFDTDBasicCoff()->getSigmaY(), xt::range(is, ie + 1),
             xt::range(js, je), ks) = sigma;
    return;
  }
}
}  // namespace xfdtd