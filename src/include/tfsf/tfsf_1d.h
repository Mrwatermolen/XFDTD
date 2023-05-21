#ifndef _TFSF_1D_H_
#define _TFSF_1D_H_

#include "tfsf/tfsf.h"

namespace xfdtd {

/**
 * @brief 1D TF/SF do not have any practical importance, but the nimplementation
 * of it is helpful in understanding more complex 2D and 3D TF/SF.
 *
 */
class TFSF1D : public TFSF {
 public:
  TFSF1D(SpatialIndex distance_z, double theta_inc, double e_i_0,
         std::unique_ptr<Waveform> waveform);
  TFSF1D(TFSF1D &&) = default;
  ~TFSF1D() override = default;

  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;

 private:
  double _exi_zn;
  double _hyi_zn;
  double _exi_zp;
  double _hyi_zp;

  double _k_dot_r0_ex_zn;
  double _k_dot_r0_hy_zn;
  double _k_dot_r0_ex_zp;
  double _k_dot_r0_hy_zp;

  void allocateKDotR() override;
  void allocateEiHi() override;
  void calculateKDotR() override;
  void caculateKDotRZN();
  void caculateKDotRZP();
};
}  // namespace xfdtd

#endif  // _TFSF_1D_H_