#ifndef _TFSF_2D_H_
#define _TFSF_2D_H_

#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {
class TFSF2D : public TFSF {
 public:
  TFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
         double ez_i0, std::unique_ptr<Waveform> waveform);
  TFSF2D(const TFSF2D &) = delete;
  TFSF2D(TFSF2D &&) = default;
  TFSF2D &operator=(const TFSF2D &) = delete;
  TFSF2D &operator=(TFSF2D &&) = delete;
  ~TFSF2D() override = default;

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

#endif  // _TFSF_2D_H_
