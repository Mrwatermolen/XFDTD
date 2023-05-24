#include "electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {
void EMF::allocateEx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ex.resize(nx, ny, nz);
  _ex.setConstant(default_value);
}

void EMF::allocateEy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ey.resize(nx, ny, nz);
  _ey.setConstant(default_value);
}

void EMF::allocateEz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _ez.resize(nx, ny, nz);
  _ez.setConstant(default_value);
}

void EMF::allocateHx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hx.resize(nx, ny, nz);
  _hx.setConstant(default_value);
}

void EMF::allocateHy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hy.resize(nx, ny, nz);
  _hy.setConstant(default_value);
}

void EMF::allocateHz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                     double default_value) {
  _hz.resize(nx, ny, nz);
  _hz.setConstant(default_value);
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