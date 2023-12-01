#include "flux/flux.h"

#include <cassert>
#include <utility>

#include "electromagnetic_field/electromagnetic_field.h"
#include "grid/grid_box.h"

namespace xfdtd {

Flux::Flux(std::unique_ptr<Shape> shape, EMComponent component)
    : _shape{std::move(shape)}, _component{component} {}

void Flux::init(std::shared_ptr<const GridSpace> grid_space,
                std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _emf = std::move(emf);
  _grid_box = _grid_space->getGridBox(_shape.get());
  if (_grid_box->getType() == GridBox::Type::CUBE ||
      _grid_box->getType() == GridBox::Type::UNDEFINED_TYPE) {
    throw std::runtime_error("The dimension of the grid box is not 1 or 2");
  }

  _value.resize({_fdtd_basic_coff->getTotalTimeStep(), _grid_box->getGridNumX(),
                 _grid_box->getGridNumY(), _grid_box->getGridNumZ()});
}

void Flux::update() {
  const auto nt{_fdtd_basic_coff->getCurrentTimeStep()};
  const auto li{_grid_box->getGridStartIndexX()};
  const auto lj{_grid_box->getGridStartIndexY()};
  const auto lk{_grid_box->getGridStartIndexZ()};
  const auto ri{_grid_box->getGridEndIndexX()};
  const auto rj{_grid_box->getGridEndIndexY()};
  const auto rk{_grid_box->getGridEndIndexZ()};

  auto dimension{_grid_box->getType()};
  if (dimension == GridBox::Type::POINT) {
  }
}

void output() {}

void Flux::captureX() {
  const auto nt{_fdtd_basic_coff->getCurrentTimeStep()};
  const auto li{_grid_box->getGridStartIndexX()};
  const auto lj{_grid_box->getGridStartIndexY()};
  const auto lk{_grid_box->getGridStartIndexZ()};
  const auto ri{_grid_box->getGridEndIndexX()};
  const auto rj{_grid_box->getGridEndIndexY()};
  const auto rk{_grid_box->getGridEndIndexZ()};
  assert(li == ri);

  for (auto j = lj; j < rj; ++j) {
    for (auto k = lk; k < rk; ++k) {
    }
  }
}

void Flux::captureY() {
  const auto nt{_fdtd_basic_coff->getCurrentTimeStep()};
  const auto li{_grid_box->getGridStartIndexX()};
  const auto lj{_grid_box->getGridStartIndexY()};
  const auto lk{_grid_box->getGridStartIndexZ()};
  const auto ri{_grid_box->getGridEndIndexX()};
  const auto rj{_grid_box->getGridEndIndexY()};
  const auto rk{_grid_box->getGridEndIndexZ()};
  assert(lj == rj);

  for (auto i = li; i < ri; ++i) {
    for (auto k = lk; k < rk; ++k) {
    }
  }
}

void Flux::captureZ() {
  const auto nt{_fdtd_basic_coff->getCurrentTimeStep()};
  const auto li{_grid_box->getGridStartIndexX()};
  const auto lj{_grid_box->getGridStartIndexY()};
  const auto lk{_grid_box->getGridStartIndexZ()};
  const auto ri{_grid_box->getGridEndIndexX()};
  const auto rj{_grid_box->getGridEndIndexY()};
  const auto rk{_grid_box->getGridEndIndexZ()};
  assert(lk == rk);

  for (auto i = li; i < ri; ++i) {
    for (auto j = lj; j < rj; ++j) {
    }
  }
}

}  // namespace xfdtd