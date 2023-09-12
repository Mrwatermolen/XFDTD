#include "boundary/boundary.h"

#include <utility>

namespace xfdtd {
void Boundary::defaultInit(std::shared_ptr<EMF> emf) {
  this->_emf = std::move(emf);
}
}  // namespace xfdtd