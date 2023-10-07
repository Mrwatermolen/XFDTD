#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
#include "shape/shape.h"
#include "util/type_define.h"

namespace xfdtd {

class Boundary {
 public:
  Boundary() = default;
  Boundary(const Boundary &) = delete;
  Boundary(Boundary &&) noexcept = default;
  Boundary &operator=(const Boundary &) = delete;
  Boundary &operator=(Boundary &&) noexcept = default;
  virtual ~Boundary() = default;

  /**
   * @brief Get the thickness of the PML
   *
   * @return SpatialIndex
   */
  virtual SpatialIndex getSize() const = 0;

  /**
   * @brief Get the normal vector. (perhaps it is not necessary)
   *
   * @return Orientation
   */
  virtual Orientation getOrientation() const = 0;

  virtual void init(std::shared_ptr<EMF> emf,
                    std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<GridSpace> grid_space,
                    std::unique_ptr<Shape> shape) = 0;

  /**
   * @brief update the H field in the boundary
   *
   */
  virtual void updateH() = 0;

  /**
   * @brief update the E field in the boundary
   *
   */
  virtual void updateE() = 0;

 protected:
  void defaultInit(std::shared_ptr<EMF> emf,
                   std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
                   std::shared_ptr<GridSpace> _grid_space,
                   std::unique_ptr<Shape> _shape);

  EMF *getEMFInstance() const;

  FDTDBasicCoff *getFDTDBasicCoff() const;

  GridSpace *getGridSpace() const;

  Shape *getShape() const;

 private:
  std::shared_ptr<EMF> _emf;
  std::shared_ptr<FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<GridSpace> _grid_space;
  std::unique_ptr<Shape> _shape;
};
}  // namespace xfdtd

#endif  // _BOUNDARY_H_
