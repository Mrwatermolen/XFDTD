#include "grid/grid_box.h"

namespace xfdtd {

GridBox::GridBox(size_t origin_i, size_t origin_j, size_t origin_k, size_t nx,
                 size_t ny, size_t nz)
    : _origin_i{origin_i},
      _origin_j{origin_j},
      _origin_k{origin_k},
      _nx{nx},
      _ny{ny},
      _nz{nz} {}

size_t GridBox::getGridNumX() const { return _nx; }

size_t GridBox::getGridNumY() const { return _ny; }

size_t GridBox::getGridNumZ() const { return _nz; }

size_t GridBox::getGridNum() const { return _nx * _ny * _nz; }

size_t GridBox::getGridOriginIndexX() const { return _origin_i; }

size_t GridBox::getGridOriginIndexY() const { return _origin_j; }

size_t GridBox::getGridOriginIndexZ() const { return _origin_k; }

size_t GridBox::getGridStartIndexX() const { return _origin_i; }

size_t GridBox::getGridStartIndexY() const { return _origin_j; }

size_t GridBox::getGridStartIndexZ() const { return _origin_k; }

size_t GridBox::getGridEndIndexX() const { return _origin_i + _nx; }

size_t GridBox::getGridEndIndexY() const { return _origin_j + _ny; }

size_t GridBox::getGridEndIndexZ() const { return _origin_k + _nz; }

size_t GridBox::getGridCenterIndexX() const { return _origin_i + _nx / 2; }

size_t GridBox::getGridCenterIndexY() const { return _origin_j + _ny / 2; }

size_t GridBox::getGridCenterIndexZ() const { return _origin_k + _nz / 2; }

}  // namespace xfdtd