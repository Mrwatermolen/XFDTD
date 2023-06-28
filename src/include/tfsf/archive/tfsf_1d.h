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
  TFSF1D(TFSF1D&&) = default;
  ~TFSF1D() override = default;

  void init(const Cube* simulation_box, double dx, double dy, double dz,
            double dt, std::unique_ptr<GridBox> tfsf_grid_box) override;
  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;
  inline double getL0() { return _l_0; }
  inline double getExi0() { return _ex_i0; }
  inline double getHyi0() { return _hy_i0; }

 private:
  double _l_0;
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

  void allocateKDotR();
  void allocateEiHi();
  void calculateKDotR();
  void caculateKDotRZN();
  void caculateKDotRZP();
};
}  // namespace xfdtd

#endif  // _TFSF_1D_H_