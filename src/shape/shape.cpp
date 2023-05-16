#include "shape/shape.h"

namespace xfdtd {
Shape::operator std::string() const { return {"Shape: "}; }
}  // namespace xfdtd
