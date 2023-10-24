#ifndef _XFDTD_GRID_SPACE_H_
#define _XFDTD_GRID_SPACE_H_

#include <cstddef>
#include <memory>
#include <xtensor.hpp>

#include "grid/grid.h"
#include "grid/grid_box.h"
#include "grid/subregion.h"
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

  double getGridSizeMinX() const;

  double getGridSizeMinY() const;

  double getGridSizeMinZ() const;

  double getGridBaseSizeX() const;

  double getGridBaseSizeY() const;

  double getGridBaseSizeZ() const;

  // double getGridSizeX() const;

  // double getGridSizeY() const;

  // double getGridSizeZ() const;

  const xt::xarray<double>& getGridSizeArrayX() const;

  const xt::xarray<double>& getGridSizeArrayY() const;

  const xt::xarray<double>& getGridSizeArrayZ() const;

  const xt::xarray<double>& getGridSizeArrayEX() const;

  const xt::xarray<double>& getGridSizeArrayEY() const;

  const xt::xarray<double>& getGridSizeArrayEZ() const;

  const xt::xarray<double>& getGridSizeArrayHX() const;

  const xt::xarray<double>& getGridSizeArrayHY() const;

  const xt::xarray<double>& getGridSizeArrayHZ() const;

  double getGridSizeX(size_t i) const;

  double getGridSizeY(size_t j) const;

  double getGridSizeZ(size_t k) const;

  double getGridSizeHX(size_t i) const;

  double getGridSizeHY(size_t j) const;

  double getGridSizeHZ(size_t k) const;

  double getGridSizeEX(size_t i) const;

  double getGridSizeEY(size_t j) const;

  double getGridSizeEZ(size_t k) const;

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

  auto getGridView(const xt::xarray<bool>& mask) const;

  auto getGridView(const Shape* shape) const;

  void calculateSpaceSize(const std::vector<std::unique_ptr<Shape>>& shapes);

  void extendBoundaryXN(double size);

  void extendBoundaryXP(double size);

  void extendBoundaryYN(double size);

  void extendBoundaryYP(double size);

  void extendBoundaryZN(double size);

  void extendBoundaryZP(double size);

  void insertSubRegion(const Subregion* subregion);

  void generateUniformGridSpace();

  bool isUniformGridSpace() const;

 private:
  // double _dx, _dy, _dz;
  double _min_x, _min_y, _min_z;
  double _max_x, _max_y, _max_z;
  double _size_x, _size_y, _size_z;
  size_t _nx{0}, _ny{0}, _nz{0};
  bool _extended_x{false}, _extended_y{false}, _extended_z{false};
  bool _1d_space, _2d_space, _3d_space;
  xt::xarray<std::shared_ptr<Grid>> _grids;
  int _material_counter{0};
  bool _is_uniform_grid_space{true};

  double _base_grid_size_x, _base_grid_size_y, _base_grid_size_z;
  double _min_grid_size_x, _min_grid_size_y, _min_grid_size_z;
  xt::xarray<double> _e_node_coord_x, _e_node_coord_y, _e_node_coord_z;
  xt::xarray<double> _h_node_coord_x, _h_node_coord_y, _h_node_coord_z;
  xt::xarray<double> _e_node_size_x, _e_node_size_y, _e_node_size_z;
  xt::xarray<double> _h_node_size_x, _h_node_size_y, _h_node_size_z;

  size_t handleTransformX(double x) const;

  size_t handleTransformY(double y) const;

  size_t handleTransformZ(double z) const;

  void resizeSpaceSize(double new_min_x, double new_min_y, double new_min_z,
                       double new_max_x, double new_max_y, double new_max_z);

  void calculateSpaceSize(Cube* cube);

  void generateGridSpace3D();

  void generateGridSpace2D();

  void generateGridSpace1D();

  void correctGridSpaceWithSubRegion(const Subregion* subregion, double base_dl,
                                     xt::xarray<double>& e_node_coord,
                                     xt::xarray<double>& h_node_coord,
                                     xt::xarray<double>& e_node_size,
                                     xt::xarray<double>& h_node_size,
                                     size_t& n);

  void correctGridSpaceCoordinate(const xt::xarray<double>& e_node_coord,
                                  xt::xarray<double>& h_node_coord,
                                  double base_dl,
                                  xt::xarray<double>& e_node_size,
                                  xt::xarray<double>& h_node_size, size_t& n);

  auto calculateTransitionRegion(double transition_start_position,
                                 double subgrid_edge_position, double base_dl,
                                 double subgrid_dl,
                                 const xt::xarray<double>& node_coord)
      -> std::tuple<size_t, double, double, xt::xarray<double>>;

  xt::xarray<double> calculateHNodeCoord(
      const xt::xarray<double>& e_node_coord);

  xt::xarray<double> calculateGridSizeEArray(
      const xt::xarray<double>& e_node_coord) const;

  xt::xarray<double> calculateGridSizeHArray(
      const xt::xarray<double>& _h_node_coord, double base_dl) const;

  /**
   * @brief
   *
   * @param node_coord
   * @param start_node
   * @param end_node
   * @param insert_node_size_arr
   */
  void insertSubRegionCoord(xt::xarray<double>& node_coord, size_t start_node,
                            size_t end_node,
                            const xt::xarray<double>& insert_node_size_arr);
};

inline auto GridSpace::getGridView(const xt::xarray<bool>& mask) const {
  return xt::filter(_grids, mask);
}

inline auto GridSpace::getGridView(const Shape* shape) const {
  return getGridView(getShapeMask(shape));
}

}  // namespace xfdtd

#endif  // _XFDTD_GRID_SPACE_H_
