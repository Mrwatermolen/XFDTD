#ifndef _TFSF_2D_H_
#define _TFSF_2D_H_

#include "tfsf/tfsf.h"
namespace xfdtd {
class TFSF2D : public TFSF {
 public:
  TFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
         double ez_i0, std::unique_ptr<Waveform> waveform);
  //   TFSF2D(const TFSF2D &) = delete;
  //   TFSF2D &operator=(const TFSF2D &) = delete;
  //   TFSF2D(TFSF2D &&) noexcept = default;
  //   TFSF2D &operator=(TFSF2D &&) = delete;
  ~TFSF2D() override = default;

  void init(const GridSpace* grid_space, const FDTDBasicCoff* fdtd_basic_coff,
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

  xt::xarray<double> _hx_inc;
  xt::xarray<double> _hy_inc;
  xt::xarray<double> _ez_inc;

  xt::xarray<double> _projection_x_full;
  xt::xarray<double> _projection_y_full;
  xt::xarray<double> _projection_x_half;
  xt::xarray<double> _projection_y_half;

  double _a{0}, _b{0};

  double getIncidentHx(int i, int j, int k);
  double getIncidentHy(int i, int j, int k);
  double getIncidentEz(int i, int j, int k);
};
}  // namespace xfdtd

#endif  // _TFSF_2D_H_