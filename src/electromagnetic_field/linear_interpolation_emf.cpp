#ifndef _XFDTD_LINEAR_INTERPOLATION_EMF_H_
#define _XFDTD_LINEAR_INTERPOLATION_EMF_H_

#include "electromagnetic_field/electromagnetic_field.h"
#include "util/type_define.h"

namespace xfdtd {

double EMF::getEMComponentAxisCenter(EMComponent c, Axis axis, SpatialIndex i,
                                     SpatialIndex j, SpatialIndex k) const {
  switch (axis) {
    case Axis::X:
      switch (c) {
        case EMComponent::EX:
          return getExAxisXCenter(i, j, k);
        case EMComponent::EY:
          return getEyAxisXCenter(i, j, k);
        case EMComponent::EZ:
          return getEzAxisXCenter(i, j, k);
        case EMComponent::HX:
          return getHxAxisXCenter(i, j, k);
        case EMComponent::HY:
          return getHyAxisXCenter(i, j, k);
        case EMComponent::HZ:
          return getHzAxisXCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Axis::Y:
      switch (c) {
        case EMComponent::EX:
          return getExAxisYCenter(i, j, k);
        case EMComponent::EY:
          return getEyAxisYCenter(i, j, k);
        case EMComponent::EZ:
          return getEzAxisYCenter(i, j, k);
        case EMComponent::HX:
          return getHxAxisYCenter(i, j, k);
        case EMComponent::HY:
          return getHyAxisYCenter(i, j, k);
        case EMComponent::HZ:
          return getHzAxisYCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Axis::Z:
      switch (c) {
        case EMComponent::EX:
          return getExAxisZCenter(i, j, k);
        case EMComponent::EY:
          return getEyAxisZCenter(i, j, k);
        case EMComponent::EZ:
          return getEzAxisZCenter(i, j, k);
        case EMComponent::HX:
          return getHxAxisZCenter(i, j, k);
        case EMComponent::HY:
          return getHyAxisZCenter(i, j, k);
        case EMComponent::HZ:
          return getHzAxisZCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    default:
      throw std::runtime_error("Invalid axis");
  }
}

double EMF::getEMComponentFaceCenter(EMComponent c, Orientation orientation,
                                     SpatialIndex i, SpatialIndex j,
                                     SpatialIndex k) const {
  switch (orientation) {
    case Orientation::XN:
      switch (c) {
        case EMComponent::EX:
          return getExFaceXCenter(i, j, k);
        case EMComponent::EY:
          return getEyFaceXCenter(i, j, k);
        case EMComponent::EZ:
          return getEzFaceXCenter(i, j, k);
        case EMComponent::HX:
          return getHxFaceXCenter(i, j, k);
        case EMComponent::HY:
          return getHyFaceXCenter(i, j, k);
        case EMComponent::HZ:
          return getHzFaceXCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Orientation::YN:
      switch (c) {
        case EMComponent::EX:
          return getExFaceYCenter(i, j, k);
        case EMComponent::EY:
          return getEyFaceYCenter(i, j, k);
        case EMComponent::EZ:
          return getEzFaceYCenter(i, j, k);
        case EMComponent::HX:
          return getHxFaceYCenter(i, j, k);
        case EMComponent::HY:
          return getHyFaceYCenter(i, j, k);
        case EMComponent::HZ:
          return getHzFaceYCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Orientation::ZN:
      switch (c) {
        case EMComponent::EX:
          return getExFaceZCenter(i, j, k);
        case EMComponent::EY:
          return getEyFaceZCenter(i, j, k);
        case EMComponent::EZ:
          return getEzFaceZCenter(i, j, k);
        case EMComponent::HX:
          return getHxFaceZCenter(i, j, k);
        case EMComponent::HY:
          return getHyFaceZCenter(i, j, k);
        case EMComponent::HZ:
          return getHzFaceZCenter(i, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Orientation::XP:
      switch (c) {
        case EMComponent::EX:
          return getExFaceXCenter(i + 1, j, k);
        case EMComponent::EY:
          return getEyFaceXCenter(i + 1, j, k);
        case EMComponent::EZ:
          return getEzFaceXCenter(i + 1, j, k);
        case EMComponent::HX:
          return getHxFaceXCenter(i + 1, j, k);
        case EMComponent::HY:
          return getHyFaceXCenter(i + 1, j, k);
        case EMComponent::HZ:
          return getHzFaceXCenter(i + 1, j, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Orientation::YP:
      switch (c) {
        case EMComponent::EX:
          return getExFaceYCenter(i, j + 1, k);
        case EMComponent::EY:
          return getEyFaceYCenter(i, j + 1, k);
        case EMComponent::EZ:
          return getEzFaceYCenter(i, j + 1, k);
        case EMComponent::HX:
          return getHxFaceYCenter(i, j + 1, k);
        case EMComponent::HY:
          return getHyFaceYCenter(i, j + 1, k);
        case EMComponent::HZ:
          return getHzFaceYCenter(i, j + 1, k);
        default:
          throw std::runtime_error("Invalid Component");
      }
    case Orientation::ZP:
      switch (c) {
        case EMComponent::EX:
          return getExFaceZCenter(i, j, k + 1);
        case EMComponent::EY:
          return getEyFaceZCenter(i, j, k + 1);
        case EMComponent::EZ:
          return getEzFaceZCenter(i, j, k + 1);
        case EMComponent::HX:
          return getHxFaceZCenter(i, j, k + 1);
        case EMComponent::HY:
          return getHyFaceZCenter(i, j, k + 1);
        case EMComponent::HZ:
          return getHzFaceZCenter(i, j, k + 1);
        default:
          throw std::runtime_error("Invalid Component");
      }
    default:
      throw std::runtime_error("Invalid orientation");
  }
}

double EMF::getEMComponentGridCenter(EMComponent c, SpatialIndex i,
                                     SpatialIndex j, SpatialIndex k) const {
  switch (c) {
    case EMComponent::EX:
      return getExGridCenter(i, j, k);
    case EMComponent::EY:
      return getEyGridCenter(i, j, k);
    case EMComponent::EZ:
      return getEzGridCenter(i, j, k);
    case EMComponent::HX:
      return getHxGridCenter(i, j, k);
    case EMComponent::HY:
      return getHyGridCenter(i, j, k);
    case EMComponent::HZ:
      return getHzGridCenter(i, j, k);
    default:
      throw std::runtime_error("Invalid Component");
  }
}

double EMF::getExAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getEx(i, j, k);
}

double EMF::getEyAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEyFaceZCenter(i, j, k) + getEyFaceZCenter(i, j - 1, k));
}

double EMF::getEzAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEzFaceYCenter(i, j, k) + getEzFaceYCenter(i, j, k - 1));
}

double EMF::getHxAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHxFaceZCenter(i, j, k) + getHxFaceZCenter(i, j - 1, k));
}

double EMF::getHyAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHy(i, j, k) + getHy(i, j, k - 1));
}

double EMF::getHzAxisXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHz(i, j, k) + getHz(i, j - 1, k));
}

double EMF::getExAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getExFaceZCenter(i, j, k) + getExFaceZCenter(i - 1, j, k));
}

double EMF::getEyAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getEy(i, j, k);
}

double EMF::getEzAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEzFaceXCenter(i, j, k) + getEzFaceXCenter(i, j, k - 1));
}

double EMF::getHxAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHx(i, j, k) + getHx(i, j, k - 1));
}

double EMF::getHyAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHyFaceXCenter(i, j, k) + getHyFaceXCenter(i, j, k - 1));
}

double EMF::getHzAxisYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHz(i, j, k) + getHz(i - 1, j, k));
}

double EMF::getExAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getExFaceYCenter(i, j, k) + getExFaceYCenter(i - 1, j, k));
}

double EMF::getEyAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEyFaceXCenter(i, j, k) + getEyFaceXCenter(i, j - 1, k));
}

double EMF::getEzAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getEz(i, j, k);
}

double EMF::getHxAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHx(i, j, k) + getHx(i, j - 1, k));
}

double EMF::getHyAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHy(i, j, k) + getHy(i - 1, j, k));
}

double EMF::getHzAxisZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getHzFaceXCenter(i, j, k) + getHzFaceXCenter(i, j - 1, k));
}

double EMF::getExFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getExGridCenter(i, j, k) + getExGridCenter(i - 1, j, k));
}

double EMF::getEyFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEy(i, j, k) + getEy(i, j, k + 1));
}

double EMF::getEzFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEz(i, j, k) + getEz(i, j + 1, k));
}

double EMF::getHxFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getHx(i, j, k);
}

// NOTE: it's ok
double EMF::getHyFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHyAxisZCenter(i, j, k) + getHyAxisZCenter(i, j + 1, k));
  return 0.25 * (getHy(i, j, k) + getHy(i, j + 1, k) + getHy(i - 1, j, k) +
                 getHy(i - 1, j + 1, k));
}

// NOTE: it's ok
double EMF::getHzFaceXCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHzAxisYCenter(i, j, k) + getHzAxisYCenter(i, j, k + 1));
  return 0.25 * (getHz(i, j, k) + getHz(i, j, k + 1) + getHz(i - 1, j, k) +
                 getHz(i - 1, j, k + 1));
}

double EMF::getExFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEx(i, j, k) + getEx(i, j + 1, k));
}

double EMF::getEyFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEyGridCenter(i, j, k) + getEyGridCenter(i, j - 1, k));
}

double EMF::getEzFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEz(i, j, k) + getEz(i + 1, j, k));
}

double EMF::getHxFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHxAxisZCenter(i, j, k) + getHxAxisZCenter(i + 1, j, k));
  return 0.25 * (getHx(i + 1, j, k) + getHx(i, j, k) + getHx(i + 1, j - 1, k) +
                 getHx(i, j - 1, k));
}

double EMF::getHyFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getHy(i, j, k);
}

double EMF::getHzFaceYCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHzAxisXCenter(i, j, k) + getHzAxisXCenter(i, j, k + 1));
  return 0.25 * (getHz(i, j, k) + getHz(i, j, k + 1) + getHz(i, j - 1, k) +
                 getHz(i, j - 1, k + 1));
}

double EMF::getExFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEx(i, j, k) + getEx(i, j, k + 1));
}

double EMF::getEyFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEy(i, j, k) + getEy(i + 1, j, k));
}

double EMF::getEzFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return 0.5 * (getEzGridCenter(i, j, k) + getEzGridCenter(i, j, k + 1));
}

double EMF::getHxFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHxAxisYCenter(i, j, k) + getHxAxisYCenter(i + 1, j, k));
  return 0.25 * (getHx(i + 1, j, k) + getHx(i, j, k) + getHx(i + 1, j, k - 1) +
                 getHx(i, j, k - 1));
}

double EMF::getHyFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  // return 0.5 * (getHyAxisXCenter(i, j, k) + getHyAxisXCenter(i, j + 1, k));
  return 0.25 * (getHy(i, j, k) + getHy(i, j + 1, k) + getHy(i, j, k - 1) +
                 getHy(i, j + 1, k - 1));
}

double EMF::getHzFaceZCenter(SpatialIndex i, SpatialIndex j,
                             SpatialIndex k) const {
  return getHz(i, j, k);
}

double EMF::getExGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.25 * (getEx(i, j, k) + getEx(i, j + 1, k) + getEx(i, j, k + 1) +
                 getEx(i, j + 1, k + 1));
}

double EMF::getEyGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.25 * (getEy(i, j, k) + getEy(i + 1, j, k) + getEy(i, j, k + 1) +
                 getEy(i + 1, j, k + 1));
}

double EMF::getEzGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.25 * (getEz(i, j, k) + getEz(i + 1, j, k) + getEz(i, j + 1, k) +
                 getEz(i + 1, j + 1, k));
}

double EMF::getHxGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.5 * (getHx(i, j, k) + getHx(i + 1, j, k));
}

double EMF::getHyGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.5 * (getHy(i, j, k) + getHy(i, j + 1, k));
}

double EMF::getHzGridCenter(SpatialIndex i, SpatialIndex j,
                            SpatialIndex k) const {
  return 0.5 * (getHz(i, j, k) + getHz(i, j, k + 1));
}

}  // namespace xfdtd

#endif  // _XFDTD_LINEAR_INTERPOLATION_EMF_H_
