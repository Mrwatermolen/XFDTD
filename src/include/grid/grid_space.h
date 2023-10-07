#ifndef _XFDTD_GRID_SPACE_H_
#define _XFDTD_GRID_SPACE_H_

#include <memory>
#include <stdexcept>
#include <xtensor.hpp>

#include "grid/grid.h"
#include "grid/grid_box.h"
#include "shape/cube.h"
#include "shape/shape.h"

namespace xfdtd {

class GridSpace {
 public:
  GridSpace(double dx, double dy, double dz);
  ~GridSpace() = default;

  size_t getGridNumX() const;

  size_t getGridNumY() const;

  size_t getGridNumZ() const;

  double getGridSizeX() const;

  double getGridSizeY() const;

  double getGridSizeZ() const;

  double getGridOriginX(size_t i) const;

  double getGridOriginY(size_t j) const;

  double getGridOriginZ(size_t k) const;

  double getGridCenterX(size_t i) const;

  double getGridCenterY(size_t j) const;

  double getGridCenterZ(size_t k) const;

  size_t getGridIndexI(double x) const;

  size_t getGridIndexJ(double y) const;

  size_t getGridIndexK(double z) const;

  Grid getGridIndex(double x, double y, double z) const;

  int getGridMaterialIndex(size_t i, size_t j, size_t k) const;

  PointVector getGridCenterPoint(size_t i, size_t j, size_t k) const;

  std::unique_ptr<Cube> getGridCube(size_t i, size_t j, size_t k) const;

  xt::xarray<bool> getShapeMask(const Shape* shape) const;

  std::unique_ptr<GridBox> getGridBox(const Shape* shape) const;

  double getGridSpaceMinX() const;

  double getGridSpaceMinY() const;

  double getGridSpaceMinZ() const;

  double getGridSpaceMaxX() const;

  double getGridSpaceMaxY() const;

  double getGridSpaceMaxZ() const;

  double getGridSpaceSizeX() const;

  double getGridSpaceSizeY() const;

  double getGridSpaceSizeZ() const;

  double getGridVolume() const;

  auto getGridView(const xt::xarray<bool>& mask) const;

  auto getGridView(const Shape* shape) const;

  void calculateSpaceSize(const std::vector<std::unique_ptr<Shape>>& shapes);

  void extendBoundaryXN(double size);

  void extendBoundaryXP(double size);

  void extendBoundaryYN(double size);

  void extendBoundaryYP(double size);

  void extendBoundaryZN(double size);

  void extendBoundaryZP(double size);

  void generateGridSpace();

 private:
  double _dx, _dy, _dz;
  double _min_x, _min_y, _min_z;
  double _max_x, _max_y, _max_z;
  double _size_x, _size_y, _size_z;
  size_t _nx{0}, _ny{0}, _nz{0};
  bool _extended_x{false}, _extended_y{false}, _extended_z{false};
  bool _1d_space, _2d_space, _3d_space;
  xt::xarray<std::shared_ptr<Grid>> _grids;
  int _material_counter{0};

  size_t handleTransformX(double x) const;

  size_t handleTransformY(double y) const;

  size_t handleTransformZ(double z) const;

  void resizeSpaceSize(double new_min_x, double new_min_y, double new_min_z,
                       double new_max_x, double new_max_y, double new_max_z);

  void calculateSpaceSize(Cube* cube);

  void generateGridSpace3D();

  void generateGridSpace2D();

  void generateGridSpace1D();
};

inline auto GridSpace::getGridView(const xt::xarray<bool>& mask) const {
  return xt::filter(_grids, mask);
}

inline auto GridSpace::getGridView(const Shape* shape) const {
  return getGridView(getShapeMask(shape));
}

}  // namespace xfdtd

#endif  // _XFDTD_GRID_SPACE_H_
