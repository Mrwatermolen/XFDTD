#ifndef _XFDTD_THIN_WIRE_H_
#define _XFDTD_THIN_WIRE_H_

#include "object/object.h"
#include "util/type_define.h"
namespace xfdtd {

class ThinWire : public Object {
 public:
  ThinWire(PointVector center, double radius, double height, Axis axis,
           std::string_view name = "");
  ~ThinWire() override = default;

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf);
  bool isPointInside(const PointVector& point) override;
  void correctCece(xt::xarray<bool>& mask, EFTA& cece, double dt);

 private:
  bool _is_init{false};
  xt::xindex _index_min{0, 0, 0};
  xt::xindex _index_max{0, 0, 0};
  void calculateIndexMin(const xt::xarray<bool>& mask);
  void calculateIndexMax(const xt::xarray<bool>& mask);
  void initThinWire(const xt::xarray<bool>& mask);
};
}  // namespace xfdtd

#endif  // _XFDTD_THIN_WIRE_H_