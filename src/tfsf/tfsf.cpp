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
      _k{Eigen::Vector3d{sin(_theta_inc) * cos(_phi_inc),
                         sin(_theta_inc) * sin(_phi_inc), cos(_theta_inc)}},
      _waveform{std::move(waveform)} {
  _ex_i0 = cos(_theta_inc) * cos(_phi_inc) * _e_theta - sin(_phi_inc) * _e_phi;
  _ey_i0 = cos(_theta_inc) * sin(_phi_inc) * _e_theta + cos(_phi_inc) * _e_phi;
  _ez_i0 = -sin(_theta_inc) * _e_theta;
  _hx_i0 =
      -(cos(_theta_inc) * cos(_phi_inc) * _e_phi + sin(_phi_inc) * _e_theta) /
      constant::ETA_0;
  _hy_i0 =
      -(cos(_theta_inc) * sin(phi_inc) * _e_phi - cos(_phi_inc) * _e_theta) /
      constant::ETA_0;
  _hz_i0 = sin(_theta_inc) * _e_phi / constant::ETA_0;
}

void TFSF::init(const Cube *simulation_box, double dx, double dy, double dz,
                double dt, TFSFBoundaryIndex tfsf_boundary_index) {
  if (simulation_box == nullptr) {
    throw std::runtime_error("simulation_box is nullptr");
  }

  initTFSF(simulation_box, dx, dy, dz, dt, tfsf_boundary_index);
  allocateKDotR();
  allocateEiHi();
  calculateKDotR();
}

void TFSF::initTFSF(const Cube *simulation_box, double dx, double dy, double dz,
                    double dt, TFSFBoundaryIndex tfsf_boundary_index) {
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

  Eigen::Matrix<double, 8, 3> fdtd_domain_matrix;
  fdtd_domain_matrix(0, 0) = simulation_box->getXmin();
  fdtd_domain_matrix(0, 1) = simulation_box->getYmin();
  fdtd_domain_matrix(0, 2) = simulation_box->getZmin();

  fdtd_domain_matrix(1, 0) = simulation_box->getXmin();
  fdtd_domain_matrix(1, 1) = simulation_box->getYmin();
  fdtd_domain_matrix(1, 2) = simulation_box->getZmax();

  fdtd_domain_matrix(2, 0) = simulation_box->getXmin();
  fdtd_domain_matrix(2, 1) = simulation_box->getYmax();
  fdtd_domain_matrix(2, 2) = simulation_box->getZmin();

  fdtd_domain_matrix(3, 0) = simulation_box->getXmin();
  fdtd_domain_matrix(3, 1) = simulation_box->getYmax();
  fdtd_domain_matrix(3, 2) = simulation_box->getZmax();

  fdtd_domain_matrix(4, 0) = simulation_box->getXmax();
  fdtd_domain_matrix(4, 1) = simulation_box->getYmin();
  fdtd_domain_matrix(4, 2) = simulation_box->getZmin();

  fdtd_domain_matrix(5, 0) = simulation_box->getXmax();
  fdtd_domain_matrix(5, 1) = simulation_box->getYmin();
  fdtd_domain_matrix(5, 2) = simulation_box->getZmax();

  fdtd_domain_matrix(6, 0) = simulation_box->getXmax();
  fdtd_domain_matrix(6, 1) = simulation_box->getYmax();
  fdtd_domain_matrix(6, 2) = simulation_box->getZmin();

  fdtd_domain_matrix(7, 0) = simulation_box->getXmax();
  fdtd_domain_matrix(7, 1) = simulation_box->getYmax();
  fdtd_domain_matrix(7, 2) = simulation_box->getZmax();

  auto k_dot_r{fdtd_domain_matrix * _k};
  auto k_dot_r_min{k_dot_r.minCoeff()};
  _l_0 = (k_dot_r_min) / constant::C_0;
}
}  // namespace xfdtd