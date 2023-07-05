#include "boundary/boundary.h"

#include <utility>

namespace xfdtd {
void Boundary::defaultInit(std::shared_ptr<EMF> _emf) {
  this->_emf = std::move(_emf);
}
}  // namespace xfdtd