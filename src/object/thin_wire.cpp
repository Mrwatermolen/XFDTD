#include "object/thin_wire.h"

#include <memory>
#include <utility>

#include "shape/cylinder.h"

namespace xfdtd {
ThinWire::ThinWire(PointVector center, double radius, double height, Axis axis,
                   std::string_view name)
    : Object(name,
             std::make_unique<Cylinder>(std::move(center), radius, height),
             std::make_unique<Material>("PEC", 1, 1, 1e10, 0)) {}

void ThinWire::correctCece(xt::xarray<bool>& mask, EFTA& cece, double dt) {
  initThinWire(mask);
  // Only for z
  auto i_min{_index_min[0]};
  auto j_min{_index_min[1]};
  auto k_min{_index_min[2]};
  auto k_max{_index_max[2]};

  xt::view(cece, i_min, j_min, xt::range(k_min, k_max)) = 0;
}

void ThinWire::initThinWire(const xt::xarray<bool>& mask) {
  if (_is_init) {
    return;
  }

  calculateIndexMin(mask);
  calculateIndexMax(mask);
  _is_init = true;
}

void ThinWire::calculateIndexMin(const xt::xarray<bool>& mask) {
  if (mask.shape().size() < 3) {
    throw std::runtime_error("mask must be 3D");
  }
  for (size_t i{0}; i < mask.shape()[0]; ++i) {
    for (size_t j{0}; j < mask.shape()[1]; ++j) {
      for (size_t k{0}; k < mask.shape()[2]; ++k) {
        if (mask(i, j, k)) {
          _index_min = {i, j, k};
          return;
        }
      }
    }
  }
}

void ThinWire::calculateIndexMax(const xt::xarray<bool>& mask) {
  if (mask.shape().size() < 3) {
    throw std::runtime_error("mask must be 3D");
  }
  for (size_t i{mask.shape()[0] - 1}; i > 0; --i) {
    for (size_t j{mask.shape()[1] - 1}; j > 0; --j) {
      for (size_t k{mask.shape()[2] - 1}; k > 0; --k) {
        if (mask(i, j, k)) {
          _index_max = {i, j, k};
          return;
        }
      }
    }
  }
}
}  // namespace xfdtd
