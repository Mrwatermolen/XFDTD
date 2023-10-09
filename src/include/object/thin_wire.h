#ifndef _XFDTD_THIN_WIRE_H_
#define _XFDTD_THIN_WIRE_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
#include "object/object.h"
#include "shape/cylinder.h"
#include "util/type_define.h"

namespace xfdtd {

/**
 * @brief PEC wire of circular cross section
 *
 */
class ThinWire : public Object {
 public:
  explicit ThinWire(const std::string& name, std::unique_ptr<Cylinder> shape);

  ~ThinWire() override = default;

  void init(int index, std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<GridSpace> grid_space,
            const std::shared_ptr<EMF>& emf) override;

  void correctFDTDCoff() override;

 private:
  Axis _axis;
  double _radius;
  void initThinWire(const xt::xarray<bool>& mask);

  bool isEnoughThin() const;
};
}  // namespace xfdtd

#endif  // _XFDTD_THIN_WIRE_H_