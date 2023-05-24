#ifndef _ELECTROMAGNETIC_FIELD_H_
#define _ELECTROMAGNETIC_FIELD_H_

#include "util/type_define.h"

namespace xfdtd {

enum class EMComponent { EX, EY, EZ, HX, HY, HZ };
class EMF {
 public:
  EMF() = default;
  inline EFTA& getEx() { return _ex; }
  inline EFTA& getEy() { return _ey; }
  inline EFTA& getEz() { return _ez; }
  inline EFTA& getHx() { return _hx; }
  inline EFTA& getHy() { return _hy; }
  inline EFTA& getHz() { return _hz; }
  inline double& getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ex(i, j, k);
  }
  inline double& getEy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ey(i, j, k);
  }
  inline double& getEz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ez(i, j, k);
  }
  inline double& getHx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _hx(i, j, k);
  }
  inline double& getHy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _hy(i, j, k);
  }
  inline double& getHz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _hz(i, j, k);
  }

  EFTA& getEMComponent(EMComponent c);
  double& getEMComponent(EMComponent c, SpatialIndex i, SpatialIndex j,
                         SpatialIndex k);

  void allocateEx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);
  void allocateEy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

  void allocateEz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

  void allocateHx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

  void allocateHy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

  void allocateHz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

 private:
  EFTA _ex;
  EFTA _ey;
  EFTA _ez;
  EFTA _hx;
  EFTA _hy;
  EFTA _hz;
};
}  // namespace xfdtd

#endif  // _ELECTROMAGNETIC_FIELD_H_