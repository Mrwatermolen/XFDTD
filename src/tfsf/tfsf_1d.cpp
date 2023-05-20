#include "tfsf/tfsf_1d.h"

#include <utility>

#include "electromagnetic.h"
#include "tfsf/tfsf.h"
#include "util/constant.h"

namespace xfdtd {
TFSF1D::TFSF1D(SpatialIndex distance_z, bool n, double e_i_0,
               std::unique_ptr<Waveform> waveform)
    : TFSF(distance_z, distance_z, distance_z,
           static_cast<int>(n) * (constant::PI), 0, e_i_0, 0,
           std::move(waveform)),
      _n{n} {
  _ex_i0 = cos(getIncidentTheta()) * getETehta();
  _hy_i0 = -(-getETehta()) / constant::ETA_0;
}

void TFSF1D::init(const Cube* simulation_box, double dx, double dy, double dz,
                  double dt, TFSF::TFSFBoundaryIndex tfsf_boundary_index) {
  if (simulation_box == nullptr) {
    throw std::runtime_error("simulation_box is nullptr");
  }

  initTFSF(simulation_box, dx, dy, dz, dt, tfsf_boundary_index);
  allocateKDotR();
  allocateEiHi();
  caculateKDotR();
}

void TFSF1D::updateIncidentField(size_t current_time_step) {
  double time_m{(current_time_step + 0.5) * getDt()};
  double time_e{(current_time_step + 1) * getDt()};

  _exi_zn =
      _ex_i0 * getIncidentFieldWaveformValueByTime(time_m - _k_dot_r0_ex_zn);
  _exi_zp =
      _ex_i0 * getIncidentFieldWaveformValueByTime(time_m - _k_dot_r0_ex_zp);
  _hyi_zn =
      _hy_i0 * getIncidentFieldWaveformValueByTime(time_e - _k_dot_r0_hy_zn);
  _hyi_zp =
      _hy_i0 * getIncidentFieldWaveformValueByTime(time_e - _k_dot_r0_hy_zp);
}

void TFSF1D::updateH() {
  if (_n) {
    auto rk{getEndIndexZ()};
    auto& hy_zp{getHy(0, 0, rk)};
    hy_zp -= (getDt() / (constant::MU_0 * getDz())) * _exi_zp;
    return;
  }

  auto lk{getStartIndexZ()};
  auto& hy_zn{getHy(0, 0, lk - 1)};
  hy_zn += (getDt() / (constant::MU_0 * getDz())) * _exi_zn;
}

void TFSF1D::updateE() {
  if (_n) {
    auto rk{getEndIndexZ()};
    auto& ex_zp{getEx(0, 0, rk)};
    ex_zp -= (getDt() / (constant::EPSILON_0 * getDz())) * _hyi_zp;
    return;
  }
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

void TFSF1D::caculateKDotR() {
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