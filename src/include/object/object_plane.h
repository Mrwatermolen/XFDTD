#ifndef _XFDTD_OBJECT_PLANE_H_
#define _XFDTD_OBJECT_PLANE_H_

#include "object/object.h"
namespace xfdtd {

/**
 * @brief class representing a plane. There must be a zero size in one dimension.
 * 
 */
class ObjectPlane : public Object {
 public:
  ObjectPlane(std::string_view name, std::unique_ptr<Shape> shape,
              std::unique_ptr<Material> material);
  ~ObjectPlane() override = default;
  void init(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<GridSpace> grid_space,
            const std::shared_ptr<EMF>& emf) override;
};
}  // namespace xfdtd

#endif  // _XFDTD_OBJECT_PLANE_H_
