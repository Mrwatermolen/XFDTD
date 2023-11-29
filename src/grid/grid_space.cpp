#include "grid/grid_space.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include "grid/grid_box.h"
#include "shape/cube.h"
#include "shape/shape.h"

namespace xfdtd {

GridSpace::GridSpace(double dx, double dy, double dz)
    : _base_grid_size_x{dx},
      _base_grid_size_y{dy},
      _base_grid_size_z{dz},
      _min_grid_size_x{dx},
      _min_grid_size_y{dy},
      _min_grid_size_z{dz} {
  _min_x = std::numeric_limits<double>::max();
  _min_y = std::numeric_limits<double>::max();
  _min_z = std::numeric_limits<double>::max();
  _max_x = std::numeric_limits<double>::min();
  _max_y = std::numeric_limits<double>::min();
  _max_z = std::numeric_limits<double>::min();
}

std::string GridSpace::toString(Dimension dimension) {
  switch (dimension) {
    case Dimension::UNDEFINED_DIMENSION:
      return "UNDEFINED_DIMENSION";
    case Dimension::ONE_DIMENSION:
      return "ONE_DIMENSION";
    case Dimension::TWO_DIMENSION:
      return "TWO_DIMENSION";
    case Dimension::THREE_DIMENSION:
      return "THREE_DIMENSION";
    default:
      return "ERROR";
  }
}

GridSpace::Dimension GridSpace::getDimension() const { return _dimension; }

size_t GridSpace::getGridNumX() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }
  if (_dimension == Dimension::TWO_DIMENSION) {
    return _nx;
  }
  if (_dimension == Dimension::THREE_DIMENSION) {
    return _nx;
  }

  throw std::runtime_error("GridSpace is not initialized");
};

size_t GridSpace::getGridNumY() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }

  if (_dimension == Dimension::TWO_DIMENSION) {
    return _ny;
  }

  if (_dimension == Dimension::THREE_DIMENSION) {
    return _ny;
  }

  throw std::runtime_error("GridSpace is not initialized");
};

size_t GridSpace::getGridNumZ() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return _nz;
  }

  if (_dimension == Dimension::TWO_DIMENSION) {
    return 1;
  }

  if (_dimension == Dimension::THREE_DIMENSION) {
    return _nz;
  }

  throw std::runtime_error("GridSpace is not initialized");
};

double GridSpace::getGridSizeMinX() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }
  return _min_grid_size_x;
}

double GridSpace::getGridSizeMinY() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }
  return _min_grid_size_y;
}

double GridSpace::getGridSizeMinZ() const {
  if (_dimension == Dimension::TWO_DIMENSION) {
    return 1;
  }

  return _min_grid_size_z;
}

double GridSpace::getGridBaseSizeX() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }

  return _base_grid_size_x;
}

double GridSpace::getGridBaseSizeY() const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 1;
  }

  return _base_grid_size_y;
}

double GridSpace::getGridBaseSizeZ() const {
  if (_dimension == Dimension::TWO_DIMENSION) {
    return 1;
  }

  return _base_grid_size_z;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayX() const {
  return _e_node_size_x;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayY() const {
  return _e_node_size_y;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayZ() const {
  return _e_node_size_z;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayEX() const {
  return _e_node_size_x;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayEY() const {
  return _e_node_size_y;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayEZ() const {
  return _e_node_size_z;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayHX() const {
  return _h_node_size_x;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayHY() const {
  return _h_node_size_y;
}

const xt::xarray<double>& GridSpace::getGridSizeArrayHZ() const {
  return _h_node_size_z;
}

double GridSpace::getGridSizeX(size_t i) const {
  // return _e_node_coord_x(i + 1) - _e_node_coord_x(i);
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 0;
  }

  return _e_node_size_x(i);
};

double GridSpace::getGridSizeY(size_t j) const {
  // return _e_node_coord_y(j + 1) - _e_node_coord_y(j);
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 0;
  }

  return _e_node_size_y(j);
};

double GridSpace::getGridSizeZ(size_t k) const {
  // return _e_node_coord_z(k + 1) - _e_node_coord_z(k);
  if (_dimension == Dimension::TWO_DIMENSION) {
    return 0;
  }

  return _e_node_size_z(k);
};

double GridSpace::getGridSizeHX(size_t i) const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return _h_node_size_x(0);
  }

  return _h_node_size_x(i);
}

double GridSpace::getGridSizeHY(size_t j) const {
  if (_dimension == Dimension::ONE_DIMENSION) {
    return 0;
  }

  return _h_node_size_y(j);
}

double GridSpace::getGridSizeHZ(size_t k) const {
  if (_dimension == Dimension::TWO_DIMENSION) {
    return 0;
  }

  return _h_node_size_z(k);
}

double GridSpace::getGridSizeEX(size_t i) const { return getGridSizeX(i); };

double GridSpace::getGridSizeEY(size_t j) const { return getGridSizeY(j); };

double GridSpace::getGridSizeEZ(size_t k) const { return getGridSizeZ(k); };

double GridSpace::getGridOriginX(size_t i) const { return _e_node_coord_x(i); };

double GridSpace::getGridOriginY(size_t j) const { return _e_node_coord_y(j); };

double GridSpace::getGridOriginZ(size_t k) const { return _e_node_coord_z(k); };

double GridSpace::getGridCenterX(size_t i) const { return _h_node_coord_x(i); };

double GridSpace::getGridCenterY(size_t j) const { return _h_node_coord_y(j); };

double GridSpace::getGridCenterZ(size_t k) const {
  if (_dimension == Dimension::TWO_DIMENSION) {
    return 0;
  }

  return _h_node_coord_z(k);
};

size_t GridSpace::getGridIndexI(double x) const { return handleTransformX(x); };

size_t GridSpace::getGridIndexJ(double y) const { return handleTransformY(y); };

size_t GridSpace::getGridIndexK(double z) const { return handleTransformZ(z); };

Grid GridSpace::getGridContainPoint(double x, double y, double z) const {
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
      PointVector{getGridSizeX(i), getGridSizeY(j), getGridSizeZ(k)});
};

xt::xarray<bool> GridSpace::getShapeMask(const Shape* shape) const {
  if (_dimension == Dimension::THREE_DIMENSION) {
    xt::xarray<bool> mask = xt::make_lambda_xfunction(
        [shape, this](const std::shared_ptr<Grid>& g) {
          return shape->isPointInside(this->getGridCenterPoint(
              g->getIndexI(), g->getIndexJ(), g->getIndexK()));
        },
        _grids);
    return mask;
  }
  if (_dimension == Dimension::TWO_DIMENSION) {
    xt::xarray<bool> mask = xt::make_lambda_xfunction(
        [shape, this](const std::shared_ptr<Grid>& g) {
          return shape->isPointInside(
              this->getGridCenterPoint(g->getIndexI(), g->getIndexJ(), 0));
        },
        _grids);
    return mask;
  }
  if (_dimension == Dimension::ONE_DIMENSION) {
    xt::xarray<bool> mask = xt::make_lambda_xfunction(
        [shape, this](const std::shared_ptr<Grid>& g) {
          return shape->isPointInside(
              this->getGridCenterPoint(g->getIndexI(), 0, 0));
        },
        _grids);
    return mask;
  }

  throw std::runtime_error("GridSpace is not initialized");
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

bool GridSpace::isUniformGridSpace() const { return _is_uniform_grid_space; }

void GridSpace::calculateSpaceSize(const Shape* shape) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }

  auto temp{shape->getWrappedBox()};
  auto box{dynamic_cast<Cube*>(temp.get())};
  if (box == nullptr) {
    throw std::runtime_error("Shape Type Error");
  }
  calculateSpaceSize(box);
}

void GridSpace::extendBoundaryXN(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_nx <= 1) {
    throw std::runtime_error("GridSpace is not extended with xn");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x - size, _min_y, _min_z, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryXP(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_nx <= 1) {
    throw std::runtime_error("GridSpace is not extended with xp");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x + size, _max_y, _max_z);
}

void GridSpace::extendBoundaryYN(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_ny <= 1) {
    throw std::runtime_error("GridSpace is not extended with yn");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x, _min_y - size, _min_z, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryYP(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_ny <= 1) {
    throw std::runtime_error("GridSpace is not extended with yp");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x, _max_y + size, _max_z);
}

void GridSpace::extendBoundaryZN(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_nz <= 1) {
    throw std::runtime_error("GridSpace is not extended with zn");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x, _min_y, _min_z - size, _max_x, _max_y, _max_z);
}

void GridSpace::extendBoundaryZP(double size) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }
  if (_dimension != Dimension::UNDEFINED_DIMENSION) {
    throw std::runtime_error("GridSpace can't be extended");
  }
  if (_nz <= 1) {
    throw std::runtime_error("GridSpace is not extended with zp");
  }
  if (size < 0) {
    throw std::runtime_error("size is less than 0");
  }

  resizeSpaceSize(_min_x, _min_y, _min_z, _max_x, _max_y, _max_z + size);
}

void GridSpace::insertSubRegion(const Subregion* subregion) {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }

  if (subregion == nullptr) {
    return;
  }
  _is_uniform_grid_space = false;

  if (subregion->getAxis() == Axis::X) {
    _min_grid_size_x = std::min(_min_grid_size_x, subregion->getDl());
    correctGridSpaceWithSubRegion(subregion, _base_grid_size_x, _e_node_coord_x,
                                  _h_node_coord_x, _e_node_size_x,
                                  _h_node_size_x, _nx);
    return;
  }

  if (subregion->getAxis() == Axis::Y) {
    _min_grid_size_y = std::min(_min_grid_size_y, subregion->getDl());
    correctGridSpaceWithSubRegion(subregion, _base_grid_size_y, _e_node_coord_y,
                                  _h_node_coord_y, _e_node_size_y,
                                  _h_node_size_y, _ny);
    return;
  }

  if (subregion->getAxis() == Axis::Z) {
    _min_grid_size_z = std::min(_min_grid_size_z, subregion->getDl());
    correctGridSpaceWithSubRegion(subregion, _base_grid_size_z, _e_node_coord_z,
                                  _h_node_coord_z, _e_node_size_z,
                                  _h_node_size_z, _nz);
    return;
  }
}

void GridSpace::generateUniformGridSpace() {
  if (_ready) {
    throw std::runtime_error("GridSpace is ready");
  }

  _nx = std::max(_nx, 1UL);
  _ny = std::max(_ny, 1UL);
  _nz = std::max(_nz, 1UL);

  auto extended_x{false};
  auto extended_y{false};
  auto extended_z{false};

  if (1 < _nx) {
    extended_x = true;
  }
  if (1 < _ny) {
    extended_y = true;
  }
  if (1 < _nz) {
    extended_z = true;
  }

  if (extended_x && extended_y && extended_z) {
    _dimension = Dimension::THREE_DIMENSION;
  } else if (extended_x && extended_x) {
    _dimension = Dimension::TWO_DIMENSION;
  } else if (extended_z) {
    _dimension = Dimension::ONE_DIMENSION;
  } else {
    _dimension = Dimension::UNDEFINED_DIMENSION;
    throw std::runtime_error("This Space can't be initialized");
  }

  generateUniformCoordinateSystem();
}

void GridSpace::tellMeOk() {
  _ready = true;

  _grids.resize({_nx, _ny, _nz});
  for (size_t i = 0; i < _nx; ++i) {
    for (size_t j = 0; j < _ny; ++j) {
      for (size_t k = 0; k < _nz; ++k) {
        _grids(i, j, k) = std::make_shared<Grid>(i, j, k);
      }
    }
  }
}

size_t GridSpace::handleTransformX(double x) const {
  if (x < _min_x) {
    throw std::out_of_range("x is out of range");
  }
  auto index{xt::argmin(xt::abs(_e_node_coord_x - x)).front()};
  return index;
  // return static_cast<size_t>(std::round<size_t>((x - _min_x) / _dx));
}

size_t GridSpace::handleTransformY(double y) const {
  if (y < _min_y) {
    throw std::out_of_range("y is out of range");
  }
  return xt::argmin(xt::abs(_e_node_coord_y - y)).front();
  // return static_cast<size_t>(std::round<size_t>((y - _min_y) / _dy));
}

size_t GridSpace::handleTransformZ(double z) const {
  if (z < _min_z) {
    throw std::out_of_range("z is out of range");
  }
  return xt::argmin(xt::abs(_e_node_coord_z - z)).front();
  // return static_cast<size_t>(std::round<size_t>((z - _min_z) / _dz));
}

void GridSpace::resizeSpaceSize(double new_min_x, double new_min_y,
                                double new_min_z, double new_max_x,
                                double new_max_y, double new_max_z) {
  if (_min_x > new_min_x) {
    _min_x = new_min_x;
  }
  if (_min_y > new_min_y) {
    _min_y = new_min_y;
  }
  if (_min_z > new_min_z) {
    _min_z = new_min_z;
  }
  if (_max_x < new_max_x) {
    _max_x = new_max_x;
  }
  if (_max_y < new_max_y) {
    _max_y = new_max_y;
  }
  if (_max_z < new_max_z) {
    _max_z = new_max_z;
  }

  size_t new_nx{
      static_cast<size_t>(std::round((_max_x - _min_x) / _base_grid_size_x))};
  size_t new_ny{
      static_cast<size_t>(std::round((_max_y - _min_y) / _base_grid_size_y))};
  size_t new_nz{
      static_cast<size_t>(std::round((_max_z - _min_z) / _base_grid_size_z))};

  if (_max_x < _min_x) {
    new_nx = 1;
  }
  if (_max_y < _min_y) {
    new_ny = 1;
  }
  if (_max_z < _min_z) {
    new_nz = 1;
  }

  _nx = std::max(_nx, new_nx);
  _ny = std::max(_ny, new_ny);
  _nz = std::max(_nz, new_nz);
}

void GridSpace::calculateSpaceSize(Cube* cube) {
  resizeSpaceSize(cube->getXmin(), cube->getYmin(), cube->getZmin(),
                  cube->getXmax(), cube->getYmax(), cube->getZmax());
}

void GridSpace::generateUniformCoordinateSystem() {
  _e_node_coord_x = _min_x + xt::arange<double>(0, _nx + 1) * _base_grid_size_x;
  _e_node_coord_y = _min_y + xt::arange<double>(0, _ny + 1) * _base_grid_size_y;
  _e_node_coord_z = _min_z + xt::arange<double>(0, _nz + 1) * _base_grid_size_z;
  correctGridSpaceCoordinateSystem(_e_node_coord_x, _h_node_coord_x,
                                   _base_grid_size_x, _e_node_size_x,
                                   _h_node_size_x, _nx);
  correctGridSpaceCoordinateSystem(_e_node_coord_y, _h_node_coord_y,
                                   _base_grid_size_y, _e_node_size_y,
                                   _h_node_size_y, _ny);
  correctGridSpaceCoordinateSystem(_e_node_coord_z, _h_node_coord_z,
                                   _base_grid_size_z, _e_node_size_z,
                                   _h_node_size_z, _nz);
}

void GridSpace::correctGridSpaceWithSubRegion(const Subregion* subregion,
                                              double base_dl,
                                              xt::xarray<double>& e_node_coord,
                                              xt::xarray<double>& h_node_coord,
                                              xt::xarray<double>& e_node_size,
                                              xt::xarray<double>& h_node_size,
                                              size_t& n) {
  auto sub_start{subregion->getStartPosition()};
  auto sub_end{subregion->getEndPosition()};
  auto sub_length{subregion->getRegionLength()};
  auto sub_transition_length{subregion->getTransitionLength()};
  auto sub_dl{subregion->getDl()};

  auto [l_tr_s_node, l_tr_s, l_tr_length, l_tr_size_arr] =
      calculateTransitionRegion(sub_start - sub_transition_length, sub_start,
                                base_dl, sub_dl, e_node_coord);
  l_tr_size_arr = xt::flip(l_tr_size_arr);

  auto [r_tr_s_node, r_tr_s, r_tr_length, r_tr_size_arr] =
      calculateTransitionRegion(sub_end + sub_transition_length, sub_end,
                                base_dl, sub_dl, e_node_coord);

  auto sub_size_arr{xt::ones<double>({sub_length / sub_dl}) * sub_dl};
  xt::xarray<double> subregion_size_arr{
      xt::concatenate(xt::xtuple(l_tr_size_arr, sub_size_arr, r_tr_size_arr))};

  // correct the grid coordinate
  insertSubRegionCoord(e_node_coord, l_tr_s_node, r_tr_s_node,
                       subregion_size_arr);

  correctGridSpaceCoordinateSystem(e_node_coord, h_node_coord, base_dl,
                                   e_node_size, h_node_size, n);
  // std::cout << "e_node_coord: " << e_node_coord << "\n";
  // std::cout << "h_node_coord: " << h_node_coord << "\n";
  // std::cout << "e_node_size: " << e_node_size << "\n";
  // std::cout << "h_node_size: " << h_node_size << "\n";
  // std::cout << "n: " << n << "\n";
}

void GridSpace::correctGridSpaceCoordinateSystem(
    const xt::xarray<double>& e_node_coord, xt::xarray<double>& h_node_coord,
    double base_dl, xt::xarray<double>& e_node_size,
    xt::xarray<double>& h_node_size, size_t& n) {
  h_node_coord = calculateHNodeCoord(e_node_coord);
  h_node_size = calculateGridSizeHArray(h_node_coord, base_dl);

  e_node_size = calculateGridSizeEArray(e_node_coord);

  n = h_node_coord.size();
}

auto GridSpace::calculateTransitionRegion(double transition_start_position,
                                          double subgrid_edge_position,
                                          double base_dl, double subgrid_dl,
                                          const xt::xarray<double>& node_coord)
    -> std::tuple<size_t, double, double, xt::xarray<double>> {
  // get closest node index
  auto closest_node_index{
      xt::argmin(xt::abs(node_coord - transition_start_position))};
  auto closest_node_position{node_coord[closest_node_index]};

  // calculate transition length
  auto transition_length{(subgrid_edge_position > closest_node_position)
                             ? (subgrid_edge_position - closest_node_position)
                             : (closest_node_position - subgrid_edge_position)};
  auto r{(transition_length + base_dl) / (transition_length + subgrid_dl)};
  auto n_transition{static_cast<size_t>(
      std::floor(std::log10(base_dl / subgrid_dl) / std::log10(r)) -
      1)};  // number of transition grids

  // calculate every grid size in transition region
  xt::xarray<double> transition_size_arr{
      subgrid_dl *
      xt::pow(r, xt::arange<double>(1, n_transition + 1))};  // can't use auto
  transition_size_arr =
      transition_size_arr * transition_length / xt::sum(transition_size_arr);

  return {closest_node_index[0], transition_start_position, transition_length,
          transition_size_arr};
}

void GridSpace::insertSubRegionCoord(
    xt::xarray<double>& node_coord, size_t start_node, size_t end_node,
    const xt::xarray<double>& insert_node_size_arr) {
  auto node_p{node_coord[start_node]};
  auto insert_node_position_arr{xt::make_lambda_xfunction(
      [&node_p](const auto& size) {
        node_p += size;
        return node_p;
      },
      insert_node_size_arr)};
  node_coord = xt::concatenate(
      xt::xtuple(xt::view(node_coord, xt::range(0, start_node + 1)),
                 insert_node_position_arr,
                 xt::view(node_coord, xt::range(end_node + 1, _))));
}

xt::xarray<double> GridSpace::calculateHNodeCoord(
    const xt::xarray<double>& e_node_coord) const {
  return (xt::view(e_node_coord, xt::range(_, -1)) +
          xt::view(e_node_coord, xt::range(1, _))) /
         2;
}

xt::xarray<double> GridSpace::calculateGridSizeEArray(
    const xt::xarray<double>& e_node_coord) const {
  return (xt::view(e_node_coord, xt::range(1, _)) -
          xt::view(e_node_coord, xt::range(_, -1)));
}

xt::xarray<double> GridSpace::calculateGridSizeHArray(
    const xt::xarray<double>& h_node_coord, double base_dl) const {
  xt::xarray<double> virtual_h_node{xt::concatenate(xt::xtuple(
      xt::xarray<double>{h_node_coord.front() - base_dl}, h_node_coord,
      xt::xarray<double>{h_node_coord.back() + base_dl}))};
  return (xt::view(virtual_h_node, xt::range(1, _)) -
          xt::view(virtual_h_node, xt::range(_, -1)));
};

}  // namespace xfdtd
