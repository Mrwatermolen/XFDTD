#ifndef _XFDTD_SUBREGION_H_
#define _XFDTD_SUBREGION_H_

#include "util/type_define.h"

namespace xfdtd {
class Subregion {
 public:
  Subregion(Axis axis, double start_position, double length,
            double transition_length, double dl);
  Axis getAxis() const;

  double getStartPosition() const;

  double getRegionLength() const;

  double getEndPosition() const;

  double getGridNum() const;

  double getTransitionLength() const;

  double getDl() const;

 private:
  Axis _axis;
  double _start_position;
  double _transition_length;
  double _dl;

  size_t _n;
  double _length, _end_position;
};
}  // namespace xfdtd

#endif  // _XFDTD_SUBREGION_H_
