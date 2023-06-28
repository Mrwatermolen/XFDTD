
#ifndef _TFSF_2D_H_
#define _TFSF_2D_H_

#include <vector>

#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {
class TFSF2D : public TFSF {
 public:
  TFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
         double ez_i0, std::unique_ptr<Waveform> waveform);
  TFSF2D(const TFSF2D &) = delete;
  TFSF2D &operator=(const TFSF2D &) = delete;
  TFSF2D(TFSF2D &&) = default;
  TFSF2D &operator=(TFSF2D &&) = delete;
  ~TFSF2D() override = default;

  void init(double dx, double dy, double dz,
            double dt, std::unique_ptr<GridBox> tfsf_grid_box) override;
  void updateIncidentField(size_t current_time_step) override;
  void updateH() override;
  void updateE() override;

  inline SpatialIndex getStrikeIndexX() const { return _strike_index_x; }
  inline SpatialIndex getStrikeIndexY() const { return _strike_index_y; }

 private:
  // fdtd update coefficient in free space
  // Ei = ceie * Ei + ceihi * Hi
  // Hi = chih * Hi + chiei * Ei
  double _ceie{1.0};
  double _ceihi;  // -(delta_t / (epsilon_0 * delta_l))
  double _chih{1.0};
  double _chiei;  // -(delta_t / (mu_0 * delta_l))

  size_t _auxiliary_array_size;

  SpatialIndex _strike_index_x;
  SpatialIndex _strike_index_y;

  // TM mode
  std::vector<double> _e_inc;
  std::vector<double> _h_inc;

  // the scalar projection on the direction k.
  std::vector<SpatialIndex> _projection_index_ez_xn;
  std::vector<SpatialIndex> _projection_index_hy_xn;

  std::vector<SpatialIndex> _projection_index_ez_yn;
  std::vector<SpatialIndex> _projection_index_hx_yn;

  std::vector<SpatialIndex> _projection_index_ez_xp;
  std::vector<SpatialIndex> _projection_index_hy_xp;

  std::vector<SpatialIndex> _projection_index_ez_yp;
  std::vector<SpatialIndex> _projection_index_hx_yp;

  // linear interpolation coefficient
  std::vector<double> _weight_ez_xn;
  std::vector<double> _weight_hy_xn;

  std::vector<double> _weight_ez_yn;
  std::vector<double> _weight_hx_yn;

  std::vector<double> _weight_ez_xp;
  std::vector<double> _weight_hy_xp;

  std::vector<double> _weight_ez_yp;
  std::vector<double> _weight_hx_yp;
};
}  // namespace xfdtd

#endif  // _TFSF_2D_H_