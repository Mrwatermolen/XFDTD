#ifndef _ELECTROMAGNETIC_FIELD_H_
#define _ELECTROMAGNETIC_FIELD_H_

#include "util/type_define.h"

namespace xfdtd {

enum class EMComponent { EX, EY, EZ, HX, HY, HZ };

/**
 * @brief A class for storing electromagnetic fields.
 *
 */
class EMF {
 public:
  EMF() = default;
  ~EMF() = default;

  const EFTA& getEx() const;
  const EFTA& getEy() const;
  const EFTA& getEz() const;
  const EFTA& getHx() const;
  const EFTA& getHy() const;
  const EFTA& getHz() const;
  const EFTA& getExPrev() const;
  const EFTA& getEyPrev() const;
  const EFTA& getEzPrev() const;

  inline EFTA& getEx() { return _ex; }
  inline EFTA& getEy() { return _ey; }
  inline EFTA& getEz() { return _ez; }
  inline EFTA& getHx() { return _hx; }
  inline EFTA& getHy() { return _hy; }
  inline EFTA& getHz() { return _hz; }
  inline EFTA& getExPrev() { return _ex_prev; }
  inline EFTA& getEyPrev() { return _ey_prev; }
  inline EFTA& getEzPrev() { return _ez_prev; }
  inline EFTA& getJx() { return _jx; }
  inline EFTA& getJy() { return _jy; }
  inline EFTA& getJz() { return _jz; }
  inline EFTA& getJxPrev() { return _jx_prev; }
  inline EFTA& getJyPrev() { return _jy_prev; }
  inline EFTA& getJzPrev() { return _jz_prev; }
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
  inline double& getExPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ex_prev(i, j, k);
  }
  inline double& getEyPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ey_prev(i, j, k);
  }
  inline double& getEzPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _ez_prev(i, j, k);
  }
  inline double& getJx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jx(i, j, k);
  }
  inline double& getJy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jy(i, j, k);
  }
  inline double& getJz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jz(i, j, k);
  }
  inline double& getJxPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jx_prev(i, j, k);
  }
  inline double& getJyPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jy_prev(i, j, k);
  }
  inline double& getJzPrev(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _jz_prev(i, j, k);
  }
  inline double getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _ex(i, j, k);
  }
  inline double getEy(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _ey(i, j, k);
  }
  inline double getEz(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _ez(i, j, k);
  }
  inline double getHx(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _hx(i, j, k);
  }
  inline double getHy(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _hy(i, j, k);
  }
  inline double getHz(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _hz(i, j, k);
  }
  inline double getExPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _ex_prev(i, j, k);
  }
  inline double getEyPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _ey_prev(i, j, k);
  }
  inline double getEzPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _ez_prev(i, j, k);
  }
  inline double getJx(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _jx(i, j, k);
  }
  inline double getJy(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _jy(i, j, k);
  }
  inline double getJz(SpatialIndex i, SpatialIndex j, SpatialIndex k) const {
    return _jz(i, j, k);
  }
  inline double getJxPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _jx_prev(i, j, k);
  }
  inline double getJyPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _jy_prev(i, j, k);
  }
  inline double getJzPrev(SpatialIndex i, SpatialIndex j,
                          SpatialIndex k) const {
    return _jz_prev(i, j, k);
  }

  const EFTA& getEMComponent(EMComponent c) const;

  EFTA& getEMComponent(EMComponent c);

  double getEMComponent(EMComponent c, SpatialIndex i, SpatialIndex j,
                        SpatialIndex k) const;

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

  void allocateExPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);
  void allocateEyPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);
  void allocateEzPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);

  void allocateJx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);
  void allocateJy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);
  void allocateJz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                  double default_value = 0.0);

  void allocateJxPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);
  void allocateJyPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);
  void allocateJzPrev(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz,
                      double default_value = 0.0);

  double getEMComponentAxisCenter(EMComponent c, Axis axis, SpatialIndex i,
                                  SpatialIndex j, SpatialIndex k) const;

  double getEMComponentFaceCenter(EMComponent c, Orientation orientation,
                                  SpatialIndex i, SpatialIndex j,
                                  SpatialIndex k) const;

  double getEMComponentGridCenter(EMComponent c, SpatialIndex i, SpatialIndex j,
                                  SpatialIndex k) const;

  double getExAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getHxAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzAxisXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getExAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getHyAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzAxisYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getExAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getHzAxisZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getExFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzFaceXCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getExFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getEyFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzFaceYCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getExFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  /**
   * @brief Be discouraged to use it
   */
  double getEzFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzFaceZCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getExGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEyGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getEzGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHxGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHyGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

  double getHzGridCenter(SpatialIndex i, SpatialIndex j, SpatialIndex k) const;

 private:
  EFTA _ex;
  EFTA _ey;
  EFTA _ez;
  EFTA _hx;
  EFTA _hy;
  EFTA _hz;
  // polarization current
  EFTA _jx, _jy, _jz;
  // polarization current arrays at previous time step
  EFTA _jx_prev, _jy_prev, _jz_prev;
  // E-field (at previous time step) for Lorentz only.
  EFTA _ex_prev, _ey_prev, _ez_prev;
};
}  // namespace xfdtd

#endif  // _ELECTROMAGNETIC_FIELD_H_