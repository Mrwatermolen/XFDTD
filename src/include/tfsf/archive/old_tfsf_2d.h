#ifndef _OLD_TFSF_2D_H_
#define _OLD_TFSF_2D_H_

#include "tfsf/archive/old_tfsf.h"
#include "util/type_define.h"

namespace xfdtd {
class OldTFSF2D : public OldTFSF {
 public:
  OldTFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
            double ez_i0, std::unique_ptr<Waveform> waveform);
  OldTFSF2D(const OldTFSF2D &) = delete;
  OldTFSF2D(OldTFSF2D &&) = default;
  OldTFSF2D &operator=(const OldTFSF2D &) = delete;
  OldTFSF2D &operator=(OldTFSF2D &&) = delete;
  ~OldTFSF2D() override = default;

  void updateIncidentField(size_t current_time_step) override;
  void updateE() override;
  void updateH() override;

 protected:
  void allocateKDotR() override;
  void allocateEiHi() override;
  void calculateKDotR() override;

 private:
  // TM mode

  DoubleArrary1D _ezi_xn;
  DoubleArrary1D _hyi_xn;

  DoubleArrary1D _ezi_yn;
  DoubleArrary1D _hxi_yn;

  DoubleArrary1D _ezi_xp;
  DoubleArrary1D _hyi_xp;

  DoubleArrary1D _ezi_yp;
  DoubleArrary1D _hxi_yp;

  DoubleArrary1D _k_dot_r0_ez_xn;
  DoubleArrary1D _k_dot_r0_hy_xn;

  DoubleArrary1D _k_dot_r0_ez_yn;
  DoubleArrary1D _k_dot_r0_hx_yn;

  DoubleArrary1D _k_dot_r0_ez_xp;
  DoubleArrary1D _k_dot_r0_hy_xp;

  DoubleArrary1D _k_dot_r0_ez_yp;
  DoubleArrary1D _k_dot_r0_hx_yp;

  void calculateKDotRXN();
  void calculateKDotRYN();
  void calculateKDotRXP();
  void calculateKDotRYP();

  void updateEi(double time);
  void updateHi(double time);

  // TE Mode
  double _hz_i0;
  double _ex_i0;
  double _ey_i0;
};
}  // namespace xfdtd

#endif  // _OLD_TFSF_2D_H_
