#include "electromagnetic.h"
#include "tfsf/archive/old_tfsf_2d.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {

OldTFSF2D::OldTFSF2D(SpatialIndex distance_x, SpatialIndex distance_y,
                     double phi_inc, double ez_i0,
                     std::unique_ptr<Waveform> waveform)
    : OldTFSF(distance_x, distance_y, 0, constant::PI / 2, phi_inc, ez_i0, 0,
              std::move(waveform)) {}

void OldTFSF2D::updateIncidentField(size_t current_time_step) {
  double time_m{(current_time_step + 0.5) * getDt()};
  double time_e{(current_time_step + 1) * getDt()};

  updateHi(time_e);
  updateEi(time_m);
}

void OldTFSF2D::updateH() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};

  for (auto i{li}; i < ri + 1; ++i) {
    auto& hx_yn{getHx(i, lj - 1, 0)};
    auto& hx_yp{getHx(i, rj, 0)};
    hx_yn += (getDt() / (constant::MU_0 * getDy())) * _ezi_yn[i - li];
    hx_yp -= (getDt() / (constant::MU_0 * getDy())) * _ezi_yp[i - li];
  }

  for (auto j{lj}; j < rj + 1; ++j) {
    auto& hy_xn{getHy(li - 1, j, 0)};
    auto hy_xp{getHy(ri, j, 0)};
    hy_xn -= (getDt() / (constant::MU_0 * getDx())) * _ezi_xn[j - lj];
    hy_xp += (getDt() / (constant::MU_0 * getDx())) * _ezi_xp[j - lj];
  }
}

void OldTFSF2D::updateE() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};

  for (auto i{li}; i < ri + 1; ++i) {
    auto& ez_yn{getEz(i, lj, 0)};
    auto& ez_yp{getEz(i, rj, 0)};
    ez_yn += (getDt() / (constant::EPSILON_0 * getDy())) * _hxi_yn[i - li];
    ez_yp -= (getDt() / (constant::EPSILON_0 * getDy())) * _hxi_yp[i - li];
  }

  for (auto j{lj}; j < rj + 1; ++j) {
    auto& ez_xn{getEz(li, j, 0)};
    auto& ez_xp{getEz(ri, j, 0)};
    ez_xn -= (getDt() / (constant::EPSILON_0 * getDx())) * _hyi_xn[j - lj];
    ez_xp += (getDt() / (constant::EPSILON_0 * getDx())) * _hyi_xp[j - lj];
  }
}

void OldTFSF2D::allocateKDotR() {
  _k_dot_r0_ez_xn = std::move(allocateDoubleArray1D(getNy() + 1));
  _k_dot_r0_hy_xn = std::move(allocateDoubleArray1D(getNy() + 1));

  _k_dot_r0_ez_yn = std::move(allocateDoubleArray1D(getNx() + 1));
  _k_dot_r0_hx_yn = std::move(allocateDoubleArray1D(getNx() + 1));

  _k_dot_r0_ez_xp = std::move(allocateDoubleArray1D(getNy() + 1));
  _k_dot_r0_hy_xp = std::move(allocateDoubleArray1D(getNy() + 1));

  _k_dot_r0_ez_yp = std::move(allocateDoubleArray1D(getNx() + 1));
  _k_dot_r0_hx_yp = std::move(allocateDoubleArray1D(getNx() + 1));
}

void OldTFSF2D::allocateEiHi() {
  _ezi_xn = std::move(allocateDoubleArray1D(getNy() + 1));
  _hyi_xn = std::move(allocateDoubleArray1D(getNy() + 1));

  _ezi_yn = std::move(allocateDoubleArray1D(getNx() + 1));
  _hxi_yn = std::move(allocateDoubleArray1D(getNx() + 1));

  _ezi_xp = std::move(allocateDoubleArray1D(getNy() + 1));
  _hyi_xp = std::move(allocateDoubleArray1D(getNy() + 1));

  _ezi_yp = std::move(allocateDoubleArray1D(getNx() + 1));
  _hxi_yp = std::move(allocateDoubleArray1D(getNx() + 1));
}

void OldTFSF2D::calculateKDotR() {
  calculateKDotRXN();
  calculateKDotRYN();
  calculateKDotRXP();
  calculateKDotRYP();
}

void OldTFSF2D::calculateKDotRXN() {
  auto k_vector{getKVector()};
  auto start_x{getTFSFCubeBox()->getPoint().x()};
  auto start_y{getTFSFCubeBox()->getPoint().y()};

  double ez_x_offset{0.0};
  double hy_x_offset{-0.5 * getDx()};
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    auto ez_y_offset{j * getDy()};
    auto hy_y_offset{j * getDy()};
    _k_dot_r0_ez_xn[j] = (k_vector.x() * (start_x + ez_x_offset) +
                          k_vector.y() * (start_y + ez_y_offset)) /
                             constant::C_0 -
                         getL0();
    _k_dot_r0_hy_xn[j] = (k_vector.x() * (start_x + hy_x_offset) +
                          k_vector.y() * (start_y + hy_y_offset)) /
                             constant::C_0 -
                         getL0();
  }
}

void OldTFSF2D::calculateKDotRXP() {
  auto k_vector{getKVector()};
  auto start_x{getTFSFCubeBox()->getPoint().x() +
               getTFSFCubeBox()->getSize().x()};
  auto start_y{getTFSFCubeBox()->getPoint().y()};

  double ez_x_offset{0.0};
  double hy_x_offset{+0.5 * getDx()};
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    auto ez_y_offset{j * getDy()};
    auto hy_y_offset{j * getDy()};
    _k_dot_r0_ez_xp[j] = (k_vector.x() * (start_x + ez_x_offset) +
                          k_vector.y() * (start_y + ez_y_offset)) /
                             constant::C_0 -
                         getL0();
    _k_dot_r0_hy_xp[j] = (k_vector.x() * (start_x + hy_x_offset) +
                          k_vector.y() * (start_y + hy_y_offset)) /
                             constant::C_0 -
                         getL0();
  }
}

void OldTFSF2D::calculateKDotRYN() {
  auto k_vector{getKVector()};
  auto start_x{getTFSFCubeBox()->getPoint().x()};
  auto start_y{getTFSFCubeBox()->getPoint().y()};

  double ez_y_offset{0.0};
  double hx_y_offset{-0.5 * getDy()};

  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    auto ez_x_offset{i * getDx()};
    auto hx_x_offset{i * getDx()};
    _k_dot_r0_ez_yn[i] = (k_vector.x() * (start_x + ez_x_offset) +
                          k_vector.y() * (start_y + ez_y_offset)) /
                             constant::C_0 -
                         getL0();
    _k_dot_r0_hx_yn[i] = (k_vector.x() * (start_x + hx_x_offset) +
                          k_vector.y() * (start_y + hx_y_offset)) /
                             constant::C_0 -
                         getL0();
  }
}

void OldTFSF2D::calculateKDotRYP() {
  auto k_vector{getKVector()};
  auto start_x{getTFSFCubeBox()->getPoint().x()};
  auto start_y{getTFSFCubeBox()->getPoint().y() +
               getTFSFCubeBox()->getSize().y()};

  auto ez_y_offset{0.0};
  auto hx_y_offset{+0.5 * getDy()};
  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    auto ez_x_offset{i * getDx()};
    auto hx_x_offset{i * getDx()};
    _k_dot_r0_ez_yp[i] = (k_vector.x() * (start_x + ez_x_offset) +
                          k_vector.y() * (start_y + ez_y_offset)) /
                             constant::C_0 -
                         getL0();
    _k_dot_r0_hx_yp[i] = (k_vector.x() * (start_x + hx_x_offset) +
                          k_vector.y() * (start_y + hx_y_offset)) /
                             constant::C_0 -
                         getL0();
  }
}

void OldTFSF2D::updateHi(double time) {
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    _hyi_xn[j] = getHyi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_hy_xn[j]);
    _hyi_xp[j] = getHyi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_hy_xp[j]);
  }

  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    _hxi_yn[i] = getHxi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_hx_yn[i]);
    _hxi_yp[i] = getHxi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_hx_yp[i]);
  }
}

void OldTFSF2D::updateEi(double time) {
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    _ezi_xn[j] = getEzi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_ez_xn[j]);
    _ezi_xp[j] = getEzi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_ez_xp[j]);
  }

  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    _ezi_yn[i] = getEzi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_ez_yn[i]);
    _ezi_yp[i] = getEzi0() *
                 getIncidentFieldWaveformValueByTime(time - _k_dot_r0_ez_yp[i]);
  }
}
}  // namespace xfdtd
