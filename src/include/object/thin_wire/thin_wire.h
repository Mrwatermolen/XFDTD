#ifndef _XFDTD_THIN_WIRE_H_
#define _XFDTD_THIN_WIRE_H_

#include "object/object.h"
#include "util/type_define.h"
namespace xfdtd {

enum class Axis { X, Y, Z };

class ThinWire : public Object {
 public:
  ThinWire(PointVector center, double radius, double height, Axis axis,
           std::string_view name = "");
  ~ThinWire() override = default;
void init();
 private:
};
}  // namespace xfdtd

#endif  // _XFDTD_THIN_WIRE_H_