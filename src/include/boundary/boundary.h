#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include <memory>
#include <utility>

#include "electromagnetic_field/electromagnetic_field.h"
#include "util/type_define.h"

namespace xfdtd {
class Simulation;

// remove it later
enum class Orientation { XN, XP, YN, YP, ZN, ZP };
class Boundary {
 public:
  Boundary() = default;
  Boundary(const Boundary &) = default;
  Boundary(Boundary &&) noexcept = default;
  Boundary &operator=(const Boundary &) = default;
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

  virtual void init(Simulation *simulation) = 0;

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
 /**
  * @brief Initialize the boundary with the EMF instance.
  * 
  * @param emf
  */
  void defaultInit(std::shared_ptr<EMF> emf);

  inline std::shared_ptr<EMF> getEMFInstance() const { return _emf; }

 private:
  std::shared_ptr<EMF> _emf;
};
}  // namespace xfdtd

#endif  // _BOUNDARY_H_
