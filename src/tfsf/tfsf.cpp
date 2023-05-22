#include "tfsf/tfsf.h"

#include <memory>

#include "shape/cube.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
TFSF::TFSF(SpatialIndex distance_x, SpatialIndex distance_y,
           SpatialIndex distance_z, double theta_inc, double phi_inc,
           double e_theta, double e_phi, std::unique_ptr<Waveform> waveform)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _theta_inc{theta_inc},
      _phi_inc{phi_inc},
      _e_theta{e_theta},
      _e_phi{e_phi},
      _k{PointVector{sin(_theta_inc) * cos(_phi_inc),
                     sin(_theta_inc) * sin(_phi_inc), cos(_theta_inc)}},
      _waveform{std::move(waveform)} {}

void TFSF::defaultInitTFSF(const Cube *simulation_box, double dx, double dy,
                           double dz, double dt,
                           TFSFBoundaryIndex tfsf_boundary_index) {
  if (simulation_box == nullptr) {
    throw std::runtime_error("TFSF::initTFSF() simulation_box is nullptr");
  }

  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;
  _tfsf_boundary_index = tfsf_boundary_index;
  _tfsf_box = std::make_unique<Cube>(
      PointVector({simulation_box->getXmin() + getStartIndexX() * getDx(),
                   simulation_box->getYmin() + getStartIndexY() * getDy(),
                   simulation_box->getZmin() + getStartIndexZ() * getDz()}),
      PointVector({getNx() * getDx(), getNy() * getDy(), getNz() * getDz()}));
}
}  // namespace xfdtd