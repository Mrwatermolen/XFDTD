#include "tfsf/tfsf_1d.h"

#include <utility>

#include "electromagnetic.h"
#include "tfsf/tfsf.h"
#include "util/constant.h"

namespace xfdtd {
TFSF1D::TFSF1D(SpatialIndex distance_z, double theta_inc, double e_i_0,
               std::unique_ptr<Waveform> waveform)
    : TFSF(0, 0, distance_z, theta_inc, 0, e_i_0, 0, std::move(waveform)) {
  _ex_i0 = cos(theta_inc) * cos(0) * e_i_0;
  _hy_i0 = -(-cos(0) * e_i_0) / constant::ETA_0;
}

void TFSF1D::init(const Cube* simulation_box, double dx, double dy, double dz,
                  double dt, TFSFBoundaryIndex tfsf_boundary_index) {
  if (simulation_box == nullptr) {
    throw std::runtime_error("TFSF::initTFSF() simulation_box is nullptr");
  }
  defaultInitTFSF(simulation_box, dx, dy, dz, dt, tfsf_boundary_index);

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

  auto k_dot_r{fdtd_domain_matrix * getKVector()};
  auto k_dot_r_min{k_dot_r.minCoeff()};
  _l_0 = (k_dot_r_min) / constant::C_0;

  allocateKDotR();
  allocateEiHi();
  calculateKDotR();
}

void TFSF1D::updateIncidentField(size_t current_time_step) {
  double time_m{(current_time_step + 0.5) * getDt()};
  double time_e{(current_time_step + 1) * getDt()};

  _exi_zn =
      getExi0() * getIncidentFieldWaveformValueByTime(time_m - _k_dot_r0_ex_zn);
  _exi_zp =
      getExi0() * getIncidentFieldWaveformValueByTime(time_m - _k_dot_r0_ex_zp);
  _hyi_zn =
      getHyi0() * getIncidentFieldWaveformValueByTime(time_e - _k_dot_r0_hy_zn);
  _hyi_zp =
      getHyi0() * getIncidentFieldWaveformValueByTime(time_e - _k_dot_r0_hy_zp);
}

void TFSF1D::updateH() {
  auto rk{getEndIndexZ()};
  auto& hy_zp{getHy(0, 0, rk)};
  hy_zp -= (getDt() / (constant::MU_0 * getDz())) * _exi_zp;

  auto lk{getStartIndexZ()};
  auto& hy_zn{getHy(0, 0, lk - 1)};
  hy_zn += (getDt() / (constant::MU_0 * getDz())) * _exi_zn;
}

void TFSF1D::updateE() {
  auto rk{getEndIndexZ()};
  auto& ex_zp{getEx(0, 0, rk)};
  ex_zp -= (getDt() / (constant::EPSILON_0 * getDz())) * _hyi_zp;

  auto lk{getStartIndexZ()};
  auto& ex_zn{getEx(0, 0, lk)};

  ex_zn += (getDt() / (constant::EPSILON_0 * getDz())) * _hyi_zn;
}

void TFSF1D::allocateKDotR() {
  _k_dot_r0_ex_zn = 0;
  _k_dot_r0_hy_zn = 0;
  _k_dot_r0_ex_zp = 0;
  _k_dot_r0_hy_zp = 0;
}

void TFSF1D::allocateEiHi() {
  _exi_zn = 0;
  _hyi_zn = 0;
  _exi_zp = 0;
  _hyi_zp = 0;
}

void TFSF1D::calculateKDotR() {
  caculateKDotRZN();
  caculateKDotRZP();
}

void TFSF1D::caculateKDotRZN() {
  auto k_vec{getKVector()};
  auto start_z{getTFSFCubeBox()->getPoint().z()};
  _k_dot_r0_ex_zn = k_vec(2) * (start_z) / constant::C_0 - getL0();
  _k_dot_r0_hy_zn =
      k_vec(2) * (start_z - 0.5 * getDz()) / constant::C_0 - getL0();
}

void TFSF1D::caculateKDotRZP() {
  auto k_vec{getKVector()};
  auto end_z{getTFSFCubeBox()->getPoint().z() +
             getTFSFCubeBox()->getSize().z()};
  _k_dot_r0_ex_zp = k_vec(2) * (end_z) / constant::C_0 - getL0();
  _k_dot_r0_hy_zp =
      k_vec(2) * (end_z + 0.5 * getDz()) / constant::C_0 - getL0();
}

}  // namespace xfdtd
