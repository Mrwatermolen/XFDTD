#ifndef _TFSF_3D_H_
#define _TFSF_3D_H_

#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {
class TFSF3D : public TFSF {
 public:
  TFSF3D(SpatialIndex distance_x, SpatialIndex distance_y,
         SpatialIndex distance_z, double e_0, double theta_inc, double phi_inc,
         double psi, std::shared_ptr<Waveform> waveform);
  TFSF3D(const TFSF3D &) = delete;
  TFSF3D &operator=(const TFSF3D &) = delete;
  TFSF3D(TFSF3D &&) = default;
  TFSF3D &operator=(TFSF3D &&) = delete;
  ~TFSF3D() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<EMF> emf) override;

  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;

 private:
  PointVector _transform_e;
  PointVector _transform_h;
  // fdtd update coefficient in free space
  // Ei = ceie * Ei + ceihi * Hi
  // Hi = chih * Hi + chiei * Ei
  double _ceie{1.0};
  double _ceihi;  // -(delta_t / (epsilon_0 * delta_l))
  double _chih{1.0};
  double _chiei;  // -(delta_t / (mu_0 * delta_l))

  size_t _auxiliary_array_size;
  double _scaled_dl;

  xt::xarray<double> _e_inc;
  xt::xarray<double> _h_inc;
  xt::xarray<double> _ex_inc;
  xt::xarray<double> _ey_inc;
  xt::xarray<double> _ez_inc;
  xt::xarray<double> _hx_inc;
  xt::xarray<double> _hy_inc;
  xt::xarray<double> _hz_inc;

  // the scalar projection on the direction k.
  xt::xarray<double> _projection_x_full;
  xt::xarray<double> _projection_y_full;
  xt::xarray<double> _projection_z_full;
  xt::xarray<double> _projection_x_half;
  xt::xarray<double> _projection_y_half;
  xt::xarray<double> _projection_z_half;

  double _a{0}, _b{0};

  double getIncidentEx(int i, int j, int k);
  double getIncidentEy(int i, int j, int k);
  double getIncidentEz(int i, int j, int k);

  double getIncidentHx(int i, int j, int k);
  double getIncidentHy(int i, int j, int k);
  double getIncidentHz(int i, int j, int k);
};
}  // namespace xfdtd

#endif  // _TFSF_3D_H_
