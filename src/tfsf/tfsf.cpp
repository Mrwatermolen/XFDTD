#include "tfsf/tfsf.h"

#include <cmath>
#include <memory>

#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
TFSF::TFSF(SpatialIndex distance_x, SpatialIndex distance_y,
           SpatialIndex distance_z, double e_0, double theta_inc,
           double phi_inc, double psi, std::unique_ptr<Waveform> waveform)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _e_0{e_0},
      _theta_inc{theta_inc},
      _phi_inc{phi_inc},
      _psi{psi},
      _sin_theta_inc{std::sin(theta_inc)},
      _cos_theta_inc{std::cos(theta_inc)},
      _sin_phi_inc{std::sin(phi_inc)},
      _cos_phi_inc{std::cos(phi_inc)},
      _sin_psi{std::sin(psi)},
      _cos_psi{std::cos(psi)},
      _k{PointVector{_sin_theta_inc * _cos_phi_inc,
                     _sin_theta_inc * _sin_phi_inc, _cos_theta_inc}},
      _waveform{std::move(waveform)} {}

void TFSF::defaultInitTFSF(double dx, double dy, double dz, double dt,
                           std::unique_ptr<GridBox> tfsf_grid_box) {
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;
  _tfsf_grid_box = std::move(tfsf_grid_box);
}
}  // namespace xfdtd