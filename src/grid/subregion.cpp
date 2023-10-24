#include "grid/subregion.h"

namespace xfdtd {
Subregion::Subregion(Axis axis, double start_position, double length,
                     double transition_length, double dl)
    : _axis{axis},
      _start_position{start_position},
      _transition_length{transition_length},
      _dl{dl} {
  _n = static_cast<size_t>(std::round(length / dl));
  _length = _n * _dl;
  _end_position = _start_position + _length;
}

Axis Subregion::getAxis() const { return _axis; }

double Subregion::getStartPosition() const { return _start_position; }

double Subregion::getRegionLength() const { return _length; }

double Subregion::getEndPosition() const { return _end_position; }

double Subregion::getGridNum() const { return _n; }

double Subregion::getTransitionLength() const { return _transition_length; }

double Subregion::getDl() const { return _dl; }

}  // namespace xfdtd
