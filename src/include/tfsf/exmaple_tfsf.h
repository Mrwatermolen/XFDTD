#ifndef _EXMAPLE_TFSF_H_
#define _EXMAPLE_TFSF_H_

#include <memory>
#include <tuple>

#include "shape/cube.h"
#include "tfsf/tfsf.h"
#include "util/type_define.h"
#include "waveform/waveform.h"

namespace xfdtd {

class ArhTFSF1D : public TFSF {
 public:
  ArhTFSF1D(SpatialIndex distance_x, SpatialIndex distance_y,
            SpatialIndex distance_z, double theta_inc, double phi_inc,
            double e_theta, double e_phi, std::unique_ptr<Waveform> waveform);

  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;

 private:
  std::unique_ptr<Waveform> _waveform;
  SpatialIndex _distance_x;
  SpatialIndex _distance_y;
  SpatialIndex _distance_z;
  double _theta_inc;
  double _phi_inc;
  double _e_theta;
  double _e_phi;
  double _ex_i0;
  double _ey_i0;
  double _ez_i0;
  double _hx_i0;
  double _hy_i0;
  double _hz_i0;

  DoubleArrary2D _eyi_xn;
  DoubleArrary2D _ezi_xn;
  DoubleArrary2D _hyi_xn;
  DoubleArrary2D _hzi_xn;
  DoubleArrary2D _ezi_yn;
  DoubleArrary2D _exi_yn;
  DoubleArrary2D _hzi_yn;
  DoubleArrary2D _hxi_yn;
  DoubleArrary2D _exi_zn;
  DoubleArrary2D _eyi_zn;
  DoubleArrary2D _hxi_zn;
  DoubleArrary2D _hyi_zn;
  DoubleArrary2D _eyi_xp;
  DoubleArrary2D _ezi_xp;
  DoubleArrary2D _hyi_xp;
  DoubleArrary2D _hzi_xp;
  DoubleArrary2D _ezi_yp;
  DoubleArrary2D _exi_yp;
  DoubleArrary2D _hzi_yp;
  DoubleArrary2D _hxi_yp;
  DoubleArrary2D _exi_zp;
  DoubleArrary2D _eyi_zp;
  DoubleArrary2D _hxi_zp;
  DoubleArrary2D _hyi_zp;

  double _dt;
  double _dx;
  double _dy;
  double _dz;
  TFSFBoundaryIndex _tfsf_boundary_index;
  Eigen::Matrix<double, 8, 3> _fdtd_domain_matrix;
  PointVector _k;
  double _l_0;

  DoubleArrary2D _k_dot_r0_ey_xn;
  DoubleArrary2D _k_dot_r0_ez_xn;
  DoubleArrary2D _k_dot_r0_hy_xn;
  DoubleArrary2D _k_dot_r0_hz_xn;

  DoubleArrary2D _k_dot_r0_ez_yn;
  DoubleArrary2D _k_dot_r0_ex_yn;
  DoubleArrary2D _k_dot_r0_hz_yn;
  DoubleArrary2D _k_dot_r0_hx_yn;

  DoubleArrary2D _k_dot_r0_ex_zn;
  DoubleArrary2D _k_dot_r0_ey_zn;
  DoubleArrary2D _k_dot_r0_hx_zn;
  DoubleArrary2D _k_dot_r0_hy_zn;

  DoubleArrary2D _k_dot_r0_ey_xp;
  DoubleArrary2D _k_dot_r0_ez_xp;
  DoubleArrary2D _k_dot_r0_hy_xp;
  DoubleArrary2D _k_dot_r0_hz_xp;

  DoubleArrary2D _k_dot_r0_ez_yp;
  DoubleArrary2D _k_dot_r0_ex_yp;
  DoubleArrary2D _k_dot_r0_hz_yp;
  DoubleArrary2D _k_dot_r0_hx_yp;

  DoubleArrary2D _k_dot_r0_ex_zp;
  DoubleArrary2D _k_dot_r0_ey_zp;
  DoubleArrary2D _k_dot_r0_hx_zp;
  DoubleArrary2D _k_dot_r0_hy_zp;

  void allocateKDotR() override;
  void allocateEiHi() override;
  void calculateKDotR() override;
  void caculateKDotRXN(double dx, double dy, double dz);
  void caculateKDotRXP(double dx, double dy, double dz);
  void caculateKDotRYN(double dx, double dy, double dz);
  void caculateKDotRYP(double dx, double dy, double dz);
  void caculateKDotRZN(double dx, double dy, double dz);
  void caculateKDotRZP(double dx, double dy, double dz);
  void updateEi(size_t current_time_step);
  void updateHi(size_t current_time_step);
  void updateEiX(double time);
  void updateEiY(double time);
  void updateEiZ(double time);
  void updateHiX(double time);
  void updateHiY(double time);
  void updateHiZ(double time);
};
}  // namespace xfdtd

#endif  // _EXMAPLE_TFSF_H_