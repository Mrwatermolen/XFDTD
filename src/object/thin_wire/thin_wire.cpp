#include "object/thin_wire/thin_wire.h"

#include <memory>
#include <utility>

#include "shape/cylinder.h"

namespace xfdtd {
ThinWire::ThinWire(PointVector center, double radius, double height, Axis axis,
                   std::string_view name)
    : Object(name,
             std::make_unique<Cylinder>(std::move(center), radius, height),
             std::make_unique<Material>("PEC", 1, 1, 1e10, 0)) {}

void ThinWire::init() {
}
}  // namespace xfdtd