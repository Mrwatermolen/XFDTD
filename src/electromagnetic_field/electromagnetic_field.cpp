#include "electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {
void EMF::allocateEx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ex = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateEy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ey = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateEz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ez = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateHx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hx = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateHy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hy = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateHz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hz = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateExPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _ex_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateEyPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _ey_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateEzPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _ez_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _jx = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _jy = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _jz = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJxPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _jx_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJyPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _jy_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

void EMF::allocateJzPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                         double default_value) {
  _jz_prev = allocateDoubleArray3D(nx, ny, nz, default_value);
}

EFTA& EMF::getEMComponent(EMComponent c) {
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

double& EMF::getEMComponent(EMComponent c, SpatialIndex i, SpatialIndex j,
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