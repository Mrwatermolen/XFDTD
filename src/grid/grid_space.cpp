#include "grid/grid_space.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <stdexcept>

#include "grid/grid_box.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/float_compare.h"

namespace xfdtd {

GridSpace::GridSpace(double dx, double dy, double dz)
    : _dx{dx}, _dy{dy}, _dz{dz} {
  _min_x = std::numeric_limits<double>::max();
  _min_y = std::numeric_limits<double>::max();
  _min_z = std::numeric_limits<double>::max();
  _max_x = std::numeric_limits<double>::min();
  _max_y = std::numeric_limits<double>::min();
  _max_z = std::numeric_limits<double>::min();
}

size_t GridSpace::getGridNumX() const { return _nx; };

size_t GridSpace::getGridNumY() const { return _ny; };

size_t GridSpace::getGridNumZ() const { return _nz; };

double GridSpace::getGridSizeX() const { return _dx; };

double GridSpace::getGridSizeY() const { return _dy; };

double GridSpace::getGridSizeZ() const { return _dz; };

double GridSpace::getGridOriginX(size_t i) const { return _min_x + (i * _dx); };

double GridSpace::getGridOriginY(size_t j) const { return _min_y + (j * _dy); };

double GridSpace::getGridOriginZ(size_t k) const { return _min_z + (k * _dz); };

double GridSpace::getGridCenterX(size_t i) const {
  return _min_x + ((i + 0.5) * _dx);
};

double GridSpace::getGridCenterY(size_t j) const {
  return _min_y + ((j + 0.5) * _dy);
};

double GridSpace::getGridCenterZ(size_t k) const {
  if (_2d_space) {
    return 0;
  }
  return _min_z + ((k + 0.5) * _dz);
};

size_t GridSpace::getGridIndexI(double x) const { return handleTransformX(x); };

size_t GridSpace::getGridIndexJ(double y) const { return handleTransformY(y); };

size_t GridSpace::getGridIndexK(double z) const { return handleTransformZ(z); };

Grid GridSpace::getGridIndex(double x, double y, double z) const {
  return {handleTransformX(x), handleTransformY(y), handleTransformZ(z)};
};

int GridSpace::getGridMaterialIndex(size_t i, size_t j, size_t k) const {
  return _grids(i, j, k)->getMaterialIndex();
}

PointVector GridSpace::getGridCenterPoint(size_t i, size_t j, size_t k) const {
  return {getGridCenterX(i), getGridCenterY(j), getGridCenterZ(k)};
};

std::unique_ptr<Cube> GridSpace::getGridCube(size_t i, size_t j,
                                             size_t k) const {
  return std::make_unique<Cube>(
      PointVector{getGridOriginX(i), getGridOriginY(j), getGridOriginZ(k)},
      PointVector{getGridSizeX(), getGridSizeY(), getGridSizeZ()});
};

xt::xarray<bool> GridSpace::getShapeMask(const Shape* shape) const {
  if (_3d_space) {
    xt::xarray<bool> mask = xt::make_lambda_xfunction(
        [shape, this](const std::shared_ptr<Grid>& g) {
          return shape->isPointInside(this->getGridCenterPoint(
              g->getIndexI(), g->getIndexJ(), g->getIndexK()));
        },
        _grids);
    return mask;
  }
  if (_2d_space) {
    xt::xarray<bool> mask = xt::make_lambda_xfunction(
        [shape, this](const std::shared_ptr<Grid>& g) {
          return shape->isPointInside(
              this->getGridCenterPoint(g->getIndexI(), g->getIndexJ(), 0));
        },
        _grids);
    return mask;
  }
  throw std::exception();
}

std::unique_ptr<GridBox> GridSpace::getGridBox(const Shape* shape) const {
  std::shared_ptr<Shape> cube_temp{shape->getWrappedBox()};
  std::shared_ptr<Cube> cube{std::dynamic_pointer_cast<Cube>(cube_temp)};
  if (cube == nullptr) {
    throw std::runtime_error("Shape Type Error");
  }
  auto start_i{handleTransformX(cube->getXmin())};
  auto start_j{handleTransformY(cube->getYmin())};
  auto start_k{handleTransformZ(cube->getZmin())};
  auto end_i{handleTransformX(cube->getXmax())};
  auto end_j{handleTransformY(cube->getYmax())};
  auto end_k{handleTransformZ(cube->getZmax())};
  if (end_i < start_i) {
    throw std::range_error("end_i < start_i");
  }
  if (end_j < start_j) {
    throw std::range_error("end_j < start_j");
  }
  if (end_k < start_k) {
    throw std::range_error("end_k < start_k");
  }
  return std::make_unique<GridBox>(start_i, start_j, start_k, end_i - start_i,
                                   end_j - start_j, end_k - start_k);
}

double GridSpace::getGridSpaceMinX() const { return _min_x; };

double GridSpace::getGridSpaceMinY() const { return _min_y; };

double GridSpace::getGridSpaceMinZ() const { return _min_z; };

double GridSpace::getGridSpaceMaxX() const { return _max_x; };

double GridSpace::getGridSpaceMaxY() const { return _max_y; };

double GridSpace::getGridSpaceMaxZ() const { return _max_z; };

double GridSpace::getGridSpaceSizeX() const { return _size_x; };

double GridSpace::getGridSpaceSizeY() const { return _size_y; };

double GridSpace::getGridSpaceSizeZ() const { return _size_z; };

double GridSpace::getGridVolume() const { return _dx * _dy * _dz; };

void GridSpace::calculateSpaceSize(
    const std::vector<std::unique_ptr<Shape>>& shapes) {
  for (auto&& e : shapes) {
    auto temp{e->getWrappedBox()};
    auto box{dynamic_cast<Cube*>(temp.get())};
    if (box == nullptr) {
      throw std::runtime_error("Shape Type Error");
    }
    calculateSpaceSize(box);
  }
}

void GridSpace::extendBoundaryXN(double size) {
  if (_nx <= 1) {
    throw std::runtime_error("GridSpace is not extended with xn");
  }
  resizeSpaceSize(_min_x - size, _min_y, _min_z, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryXP(double size) {
  if (_nx <= 1) {
    throw std::runtime_error("GridSpace is not extended with xp");
  }
  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x + size, _max_y, _max_z);
}

void GridSpace::extendBoundaryYN(double size) {
  if (_ny <= 1) {
    throw std::runtime_error("GridSpace is not extended with yn");
  }
  resizeSpaceSize(_min_x, _min_y - size, _min_z, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryYP(double size) {
  if (_ny <= 1) {
    throw std::runtime_error("GridSpace is not extended with yp");
  }
  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x, _max_y + size, _max_z);
}

void GridSpace::extendBoundaryZN(double size) {
  if (_nz <= 1) {
    throw std::runtime_error("GridSpace is not extended with zn");
  }
  resizeSpaceSize(_min_x, _min_y, _min_z - size, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryZP(double size) {
  if (_nz <= 1) {
    throw std::runtime_error("GridSpace is not extended with zp");
  }
  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x, _max_y, _max_z + size);
}

void GridSpace::generateGridSpace() {
  _nx = std::max(_nx, 1UL);
  _ny = std::max(_ny, 1UL);
  _nz = std::max(_nz, 1UL);

  if (_nx != 1) {
    _extended_x = true;
  }
  if (_ny != 1) {
    _extended_y = true;
  }
  if (_nz != 1) {
    _extended_z = true;
  }
  if (_extended_x && _extended_y && _extended_z) {
    _3d_space = true;
  } else if (_extended_x && _extended_y) {
    _2d_space = true;
  } else if (_extended_z) {
    _1d_space = true;
  } else {
    throw std::runtime_error("GridSpace is not extended");
  }
  if (_3d_space) {
    generateGridSpace3D();
  } else if (_2d_space) {
    generateGridSpace2D();
  } else if (_1d_space) {
    generateGridSpace1D();
  }
}

size_t GridSpace::handleTransformX(double x) const {
  if (x < _min_x) {
    throw std::out_of_range("x is out of range");
  }
  return static_cast<size_t>(std::round<size_t>((x - _min_x) / _dx));
}
size_t GridSpace::handleTransformY(double y) const {
  if (y < _min_y) {
    throw std::out_of_range("y is out of range");
  }
  return static_cast<size_t>(std::round<size_t>((y - _min_y) / _dy));
}
size_t GridSpace::handleTransformZ(double z) const {
  if (z < _min_z) {
    throw std::out_of_range("z is out of range");
  }
  return static_cast<size_t>(std::round<size_t>((z - _min_z) / _dz));
}

void GridSpace::resizeSpaceSize(double new_min_x, double new_min_y,
                                double new_min_z, double new_max_x,
                                double new_max_y, double new_max_z) {
  if (isGreaterOrEqual(_min_x, new_min_x)) {
    _min_x = new_min_x;
  }
  if (isGreaterOrEqual(_min_y, new_min_y)) {
    _min_y = new_min_y;
  }
  if (isGreaterOrEqual(_min_z, new_min_z)) {
    _min_z = new_min_z;
  }
  if (isLessOrEqual(_max_x, new_max_x)) {
    _max_x = new_max_x;
  }
  if (isLessOrEqual(_max_y, new_max_y)) {
    _max_y = new_max_y;
  }
  if (isLessOrEqual(_max_z, new_max_z)) {
    _max_z = new_max_z;
  }
  size_t new_nx{static_cast<size_t>(std::round((_max_x - _min_x) / _dx))};
  size_t new_ny{static_cast<size_t>(std::round((_max_y - _min_y) / _dy))};
  size_t new_nz{static_cast<size_t>(std::round((_max_z - _min_z) / _dz))};
  _nx = std::max(_nx, new_nx);
  _ny = std::max(_ny, new_ny);
  _nz = std::max(_nz, new_nz);
}

void GridSpace::calculateSpaceSize(Cube* cube) {
  resizeSpaceSize(cube->getXmin(), cube->getYmin(), cube->getZmin(),
                  cube->getXmax(), cube->getYmax(), cube->getZmax());
}

void GridSpace::generateGridSpace3D() {
  _grids.resize({_nx, _ny, _nz});
  for (size_t i = 0; i < _nx; ++i) {
    for (size_t j = 0; j < _ny; ++j) {
      for (size_t k = 0; k < _nz; ++k) {
        _grids(i, j, k) = std::make_shared<Grid>(i, j, k);
      }
    }
  }
}

void GridSpace::generateGridSpace2D() {
  _grids.resize({_nx, _ny, _nz});
  for (size_t i = 0; i < _nx; ++i) {
    for (size_t j = 0; j < _ny; ++j) {
      _grids(i, j, 0) = std::make_shared<Grid>(i, j, 0);
    }
  }
}

void GridSpace::generateGridSpace1D() {
  _grids.resize({_nx});
  for (size_t k = 0; k < _nx; ++k) {
    _grids(0, 0, k) = std::make_shared<Grid>(0, 0, 1);
  }
}

}  // namespace xfdtd