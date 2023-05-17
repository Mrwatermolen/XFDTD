#ifndef _ELECTROMAGNETIC_H_
#define _ELECTROMAGNETIC_H_

#include "util/type_define.h"

namespace xfdtd {

enum class EMComponent { EX, EY, EZ, HX, HY, HZ };

inline EFTA global_ex, global_ey, global_ez, global_hx, global_hy, global_hz;

inline void allocateEx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_ex.resize(nx, ny, nz);
  global_ex.setConstant(default_value);
}

inline void allocateEy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_ey.resize(nx, ny, nz);
  global_ey.setConstant(default_value);
}

inline void allocateEz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_ez.resize(nx, ny, nz);
  global_ez.setConstant(default_value);
}

inline void allocateHx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_hx.resize(nx, ny, nz);
  global_hx.setConstant(default_value);
}

inline void allocateHy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_hy.resize(nx, ny, nz);
  global_hy.setConstant(default_value);
}

inline void allocateHz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                       double default_value = 0.0) {
  global_hz.resize(nx, ny, nz);
  global_hz.setConstant(default_value);
}

inline DoubleArrary3D allocateDoubleArray3D(SpatialIndex nx, SpatialIndex ny,
                                            SpatialIndex nz,
                                            double default_value = 0.0) {
  auto arr{DoubleArrary3D(nx, ny, nz)};
  arr.setConstant(default_value);
  return arr;
}

inline DoubleArrary2D allocateDoubleArray2D(SpatialIndex nx, SpatialIndex ny,
                                            double default_value = 0.0) {
  auto arr{DoubleArrary2D(nx, ny)};
  arr.setConstant(default_value);
  return arr;
}

inline DoubleArrary1D allocateDoubleArray1D(SpatialIndex nz,
                                            double default_value = 0.0) {
  auto arr{DoubleArrary1D{}};
  arr.resize(nz);
  for (auto&& e : arr) {
    e = default_value;
  }
  return arr;
}

inline EFTA& getEx() { return global_ex; }
inline EFTA& getEy() { return global_ey; };
inline EFTA& getEz() { return global_ez; };
inline EFTA& getHx() { return global_hx; }
inline EFTA& getHy() { return global_hy; };
inline EFTA& getHz() { return global_hz; };

inline double getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_ex(i, j, k);
}

inline double getEy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_ey(i, j, k);
}

inline double getEz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_ez(i, j, k);
}

inline double getHx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_hx(i, j, k);
}

inline double getHy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_hy(i, j, k);
}

inline double getHz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
  return global_hz(i, j, k);
}

inline EFTA& getEMComponent(EMComponent c) {
  switch (c) {
    case EMComponent::EX:
      return getEx();
    case EMComponent::EY:
      return getEy();
    case EMComponent::EZ:
      return getEz();
    case EMComponent::HX:
      return getHx();
    case EMComponent::HY:
      return getHy();
    case EMComponent::HZ:
      return getHz();
  }
}

inline double getEMComponent(EMComponent c, SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) {
  switch (c) {
    case EMComponent::EX:
      return getEx(i, j, k);
    case EMComponent::EY:
      return getEy(i, j, k);
    case EMComponent::EZ:
      return getEz(i, j, k);
    case EMComponent::HX:
      return getHx(i, j, k);
    case EMComponent::HY:
      return getHy(i, j, k);
    case EMComponent::HZ:
      return getHz(i, j, k);
  }
}

}  // namespace xfdtd

#endif  // _ELECTROMAGNETIC_H_
