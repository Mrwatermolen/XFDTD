#ifndef _TFSF_1D_H_
#define _TFSF_1D_H_

#include "tfsf/tfsf.h"

namespace xfdtd {

class TFSF1D : public TFSF {
 public:
  TFSF1D(SpatialIndex distance_z, bool n, double e_i_0,
         std::unique_ptr<Waveform> waveform);
  TFSF1D(TFSF1D &&) = default;
  ~TFSF1D() override = default;

  void init(const Cube *simulation_box, double dx, double dy, double dz,
            double dt, TFSF::TFSFBoundaryIndex tfsf_boundary_index) override;

  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;

 private:
  bool _n;
  double _ex_i0;
  double _hy_i0;

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
  void caculateKDotR() override;
  void caculateKDotRZN();
  void caculateKDotRZP();
};
}  // namespace xfdtd

#endif  // _TFSF_1D_H_