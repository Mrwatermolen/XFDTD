#include "tfsf/exmaple_tfsf.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "electromagnetic.h"
#include "shape/cube.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
ArhTFSF1D::ArhTFSF1D(SpatialIndex distance_x, SpatialIndex distance_y,
                     SpatialIndex distance_z, double theta_inc, double phi_inc,
                     double e_theta, double e_phi,
                     std::unique_ptr<Waveform> waveform)
    : TFSF(distance_z, distance_z, distance_z, theta_inc, phi_inc, e_theta,
           e_phi, std::move(waveform)) {
  _k = Eigen::Vector3d{sin(_theta_inc) * cos(_phi_inc),
                       sin(_theta_inc) * sin(_phi_inc), cos(_theta_inc)};
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

void ArhTFSF1D::init(const Cube* simulation_box, double dx, double dy,
                     double dz, double dt,
                     TFSFBoundaryIndex tfsf_boundary_index) {
  if (simulation_box == nullptr) {
    throw std::runtime_error("simulation_box is nullptr");
  }

  allocateKDotR();
  allocateEiHi();
  caculateKDotR();
}

void ArhTFSF1D::updateIncidentField(size_t current_time_step) {
  updateEi(current_time_step);
  updateHi(current_time_step);
}

void ArhTFSF1D::updateH() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};
  auto lk{getStartIndexZ()};
  auto rk{getEndIndexZ()};

  if (getNx() == 1 && getNy() == 1) {
    // 1d
    auto& hy_zn{getHy(0, 0, lk - 1)};
    auto& hy_zp{getHy(0, 0, rk)};
    hy_zn += (_dt / (constant::MU_0 * _dz)) * _exi_zn(0, 0);
    // hy_zp -= (_dt / (constant::MU_0 * _dz)) * _exi_zp(0, 0);
    return;
  }

  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& hz_xn{getHz(li - 1, j, k)};
      auto& hz_xp{getHz(ri, j, k)};
      hz_xn += (_dt / (constant::MU_0 * _dx)) * _eyi_xn(j - lj, k - lk);
      hz_xp -= (_dt / (constant::MU_0 * _dx)) * _eyi_xp(j - lj, k - lk);
    }
  }

  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& hy_xn{getHy(li - 1, j, k)};
      auto& hy_xp{getHy(ri, j, k)};
      hy_xn -= (_dt / (constant::MU_0 * _dx)) * _ezi_xn(j - lj, k - lk);
      hy_xp += (_dt / (constant::MU_0 * _dx)) * _ezi_xp(j - lj, k - lk);
    }
  }

  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri + 1; ++i) {
      auto& hx_yn{getHx(i, lj - 1, k)};
      auto& hx_yp{getHx(i, rj, k)};
      hx_yn += (_dt / (constant::MU_0 * _dy)) * _ezi_yn(k - lk, i - li);
      hx_yp -= (_dt / (constant::MU_0 * _dy)) * _ezi_yp(k - lk, i - li);
    }
  }

  for (SpatialIndex k{lk}; k < rk + 1; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      auto& hz_yn{getHz(i, lj - 1, k)};
      auto& hz_yp{getHz(i, rj, k)};
      hz_yn -= (_dt / (constant::MU_0 * _dy)) * _exi_yn(k - lk, i - li);
      hz_yp += (_dt / (constant::MU_0 * _dy)) * _exi_yp(k - lk, i - li);
    }
  }

  for (SpatialIndex i{li}; i < ri; ++i) {
    for (SpatialIndex j{lj}; j < rj + 1; ++j) {
      auto& hy_zn{getHy(i, j, lk - 1)};
      auto& hy_zp{getHy(i, j, rk)};
      hy_zn += (_dt / (constant::MU_0 * _dz)) * _exi_zn(i - li, j - lj);
      hy_zp -= (_dt / (constant::MU_0 * _dz)) * _exi_zp(i - li, j - lj);
    }
  }

  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    for (SpatialIndex j{lj}; j < rj; ++j) {
      auto& hx_zn{getHx(i, j, lk - 1)};
      auto& hx_zp{getHx(i, j, rk)};
      hx_zn -= (_dt / (constant::MU_0 * _dz)) * _eyi_zn(i - li, j - lj);
      hx_zp += (_dt / (constant::MU_0 * _dz)) * _eyi_zp(i - li, j - lj);
    }
  }
}

void ArhTFSF1D::updateE() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};
  auto lk{getStartIndexZ()};
  auto rk{getEndIndexZ()};

  if (getNx() == 1 && getNy() == 1) {
    auto& ex_zn{getEx(0, 0, lk)};
    auto& ex_zp{getEx(0, 0, rk)};
    ex_zn += (_dt / (constant::EPSILON_0 * _dz)) * _hyi_zn(0, 0);
    // ex_zp -= (_dt / (constant::EPSILON_0 * _dz)) * _hyi_zp(0, 0);
    return;
  }

  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& ey_xn{getEy(li, j, k)};
      auto& ey_xp{getEy(ri, j, k)};
      ey_xn += (_dt / (constant::EPSILON_0 * _dx)) * _hzi_xn(j - lj, k - lk);
      ey_xp -= (_dt / (constant::EPSILON_0 * _dx)) * _hzi_xp(j - lj, k - lk);
    }
  }

  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& ez_xn{getEz(li, j, k)};
      auto& ez_xp{getEz(ri, j, k)};
      ez_xn -= (_dt / (constant::EPSILON_0 * _dx)) * _hyi_xn(j - lj, k - lk);
      ez_xp += (_dt / (constant::EPSILON_0 * _dx)) * _hyi_xp(j - lj, k - lk);
    }
  }

  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri + 1; ++i) {
      auto& ez_yn{getEz(i, lj, k)};
      auto& ez_yp{getEz(i, rj, k)};
      ez_yn += (_dt / (constant::EPSILON_0 * _dy)) * _hxi_yn(k - lk, i - li);
      ez_yp -= (_dt / (constant::EPSILON_0 * _dy)) * _hxi_yp(k - lk, i - li);
    }
  }

  for (SpatialIndex k{lk}; k < rk + 1; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      auto& ex_yn{getEx(i, lj, k)};
      auto& ex_yp{getEx(i, rj, k)};
      ex_yn -= (_dt / (constant::EPSILON_0 * _dy)) * _hzi_yn(k - lk, i - li);
      ex_yp += (_dt / (constant::EPSILON_0 * _dy)) * _hzi_yp(k - lk, i - li);
    }
  }

  for (SpatialIndex i{li}; i < ri; ++i) {
    for (SpatialIndex j{lj}; j < rj + 1; ++j) {
      auto& ex_zn{getEx(i, j, lk)};
      auto& ex_zp{getEx(i, j, rk)};
      ex_zn += (_dt / (constant::EPSILON_0 * _dz)) * _hyi_zn(i - li, j - lj);
      ex_zp -= (_dt / (constant::EPSILON_0 * _dz)) * _hyi_zp(i - li, j - lj);
    }
  }

  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    for (SpatialIndex j{lj}; j < rj; ++j) {
      auto& ey_zn{getEy(i, j, lk)};
      auto& ey_zp{getEy(i, j, rk)};
      ey_zn -= (_dt / (constant::EPSILON_0 * _dz)) * _hxi_zn(i, j);
      ey_zp += (_dt / (constant::EPSILON_0 * _dz)) * _hxi_zp(i, j);
    }
  }
}

void ArhTFSF1D::updateEi(size_t current_time_step) {
  double time{(current_time_step + 0.5) * _dt};
  updateEiX(time);
  updateEiY(time);
  updateEiZ(time);
}

void ArhTFSF1D::updateHi(size_t current_time_step) {
  double time{(current_time_step + 1) * _dt};
  updateHiX(time);
  updateHiY(time);
  updateHiZ(time);
}

void ArhTFSF1D::allocateKDotR() {
  // x
  // y z
  _k_dot_r0_ey_xn = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));
  _k_dot_r0_ez_xn = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _k_dot_r0_hy_xn = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _k_dot_r0_hz_xn = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));

  _k_dot_r0_ey_xp = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));
  _k_dot_r0_ez_xp = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _k_dot_r0_hy_xp = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _k_dot_r0_hz_xp = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));

  _k_dot_r0_ez_yn = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));
  _k_dot_r0_ex_yn = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _k_dot_r0_hz_yn = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _k_dot_r0_hx_yn = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));

  _k_dot_r0_ez_yp = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));
  _k_dot_r0_ex_yp = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _k_dot_r0_hz_yp = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _k_dot_r0_hx_yp = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));

  _k_dot_r0_ex_zn = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
  _k_dot_r0_ey_zn = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _k_dot_r0_hx_zn = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _k_dot_r0_hy_zn = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));

  _k_dot_r0_ex_zp = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
  _k_dot_r0_ey_zp = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _k_dot_r0_hx_zp = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _k_dot_r0_hy_zp = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
}

void ArhTFSF1D::allocateEiHi() {
  _eyi_xn = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));
  _ezi_xn = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _hyi_xn = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _hzi_xn = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));

  _eyi_xp = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));
  _ezi_xp = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _hyi_xp = std::move(allocateDoubleArray2D(getNy() + 1, getNz()));
  _hzi_xp = std::move(allocateDoubleArray2D(getNy(), getNz() + 1));

  _ezi_yn = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));
  _exi_yn = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _hzi_yn = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _hxi_yn = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));

  _ezi_yp = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));
  _exi_yp = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _hzi_yp = std::move(allocateDoubleArray2D(getNz() + 1, getNx()));
  _hxi_yp = std::move(allocateDoubleArray2D(getNz(), getNx() + 1));

  _exi_zn = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
  _eyi_zn = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _hxi_zn = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _hyi_zn = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));

  _exi_zp = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
  _eyi_zp = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _hxi_zp = std::move(allocateDoubleArray2D(getNx() + 1, getNy()));
  _hyi_zp = std::move(allocateDoubleArray2D(getNx(), getNy() + 1));
}

void ArhTFSF1D::caculateKDotR() {
  auto dx{getDx()};
  auto dy{getDy()};
  auto dz{getDz()};
  caculateKDotRXN(dx, dy, dz);
  caculateKDotRXP(dx, dy, dz);
  caculateKDotRYN(dx, dy, dz);
  caculateKDotRYP(dx, dy, dz);
  caculateKDotRZN(dx, dy, dz);
  caculateKDotRZP(dx, dy, dz);
}

void ArhTFSF1D::caculateKDotRXN(double dx, double dy, double dz) {
  auto start_x{_fdtd_domain_matrix(0, 0) + dx * getStartIndexX()};
  auto start_y{_fdtd_domain_matrix(0, 1) + dy * getStartIndexY()};
  auto start_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};

  double ey_x_offset{0.0};
  double ez_x_offset{0.0};
  double hy_x_offset{-0.5 * dx};
  double hz_x_offset{-0.5 * dx};

  for (SpatialIndex j{0}; j < getNy(); ++j) {
    auto ey_y_offset{(j + 0.5) * dy};
    auto hz_y_offset{(j + 0.5) * dy};
    for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
      auto ey_z_offset{k * dz};
      auto hz_z_offset{k * dz};
      _k_dot_r0_ey_xn(j, k) =
          (_k(0) * (start_x + ey_x_offset) + _k(1) * (start_y + ey_y_offset) +
           _k(2) * (start_z + ey_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hz_xn(j, k) =
          (_k(0) * (start_x + hz_x_offset) + _k(1) * (start_y + hz_y_offset) +
           _k(2) * (start_z + hz_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    auto ez_y_offset{j * dy};
    auto hy_y_offset{j * dy};
    for (SpatialIndex k{0}; k < getNz(); ++k) {
      auto ez_z_offset{(k + 0.5) * dz};
      auto hy_z_offset{(k + 0.5) * dz};
      _k_dot_r0_ez_xn(j, k) =
          (_k(0) * (start_x + ez_x_offset) + _k(1) * (start_y + ez_y_offset) +
           _k(2) * (start_z + ez_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hy_xn(j, k) =
          (_k(0) * (start_x + hy_x_offset) + _k(1) * (start_y + hy_y_offset) +
           _k(2) * (start_z + hy_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }
}

void ArhTFSF1D::caculateKDotRXP(double dx, double dy, double dz) {
  auto start_x{_fdtd_domain_matrix(0, 0) + dx * getEndIndexX()};
  auto start_y{_fdtd_domain_matrix(0, 1) + dy * getStartIndexY()};
  auto start_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};

  double ey_x_offset{0.0};
  double ez_x_offset{0.0};
  double hy_x_offset{0.5 * dx};
  double hz_x_offset{0.5 * dx};

  for (SpatialIndex j{0}; j < getNy(); ++j) {
    auto ey_y_offset{(j + 0.5) * dy};
    auto hz_y_offset{(j + 0.5) * dy};
    for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
      auto ey_z_offset{k * dz};
      auto hz_z_offset{k * dz};
      _k_dot_r0_ey_xp(j, k) =
          (_k(0) * (start_x + ey_x_offset) + _k(1) * (start_y + ey_y_offset) +
           _k(2) * (start_z + ey_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hz_xp(j, k) =
          (_k(0) * (start_x + hz_x_offset) + _k(1) * (start_y + hz_y_offset) +
           _k(2) * (start_z + hz_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    auto ez_y_offset{j * dy};
    auto hy_y_offset{j * dy};
    for (SpatialIndex k{0}; k < getNz(); ++k) {
      auto ez_z_offset{(k + 0.5) * dz};
      auto hy_z_offset{(k + 0.5) * dz};
      _k_dot_r0_ez_xp(j, k) =
          (_k(0) * (start_x + ez_x_offset) + _k(1) * (start_y + ez_y_offset) +
           _k(2) * (start_z + ez_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hy_xp(j, k) =
          (_k(0) * (start_x + hy_x_offset) + _k(1) * (start_y + hy_y_offset) +
           _k(2) * (start_z + hy_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }
}

void ArhTFSF1D::caculateKDotRYN(double dx, double dy, double dz) {
  auto min_x{_fdtd_domain_matrix(0, 0) + dx * getStartIndexX()};
  auto min_y{_fdtd_domain_matrix(0, 1) + dy * getStartIndexY()};
  auto min_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};

  double ez_y_offset{0.0};
  double ex_y_offset{0.0};
  double hz_y_offset{-0.5 * dy};
  double hx_y_offset{-0.5 * dy};

  for (SpatialIndex k{0}; k < getNz(); ++k) {
    auto ez_z_offset{(k + 0.5) * dz};
    auto hx_z_offset{(k + 0.5) * dz};
    for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
      auto ez_x_offset{i * dx};
      auto hx_x_offset{i * dx};
      _k_dot_r0_ez_yn(k, i) =
          (_k(0) * (min_x + ez_x_offset) + _k(1) * (min_y + ez_y_offset) +
           _k(2) * (min_z + ez_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hx_yn(k, i) =
          (_k(0) * (min_x + hx_x_offset) + _k(1) * (min_y + hx_y_offset) +
           _k(2) * (min_z + hx_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
    auto ex_z_offset{k * dz};
    auto hz_z_offset{k * dz};
    for (SpatialIndex i{0}; i < getNx(); ++i) {
      auto ex_x_offset{(i + 0.5) * dx};
      auto hz_x_offset{(i + 0.5) * dx};
      _k_dot_r0_ex_yn(k, i) =
          (_k(0) * (min_x + ex_x_offset) + _k(1) * (min_y + ex_y_offset) +
           _k(2) * (min_z + ex_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hz_yn(k, i) =
          (_k(0) * (min_x + hz_x_offset) + _k(1) * (min_y + hz_y_offset) +
           _k(2) * (min_z + hz_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }
}

void ArhTFSF1D::caculateKDotRYP(double dx, double dy, double dz) {
  auto max_x{_fdtd_domain_matrix(0, 0) + dx * getStartIndexX()};
  auto max_y{_fdtd_domain_matrix(0, 1) + dy * getEndIndexY()};
  auto max_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};

  double ez_y_offset{0.0};
  double ex_y_offset{0.0};
  double hz_y_offset{0.5 * dy};
  double hx_y_offset{0.5 * dy};

  for (SpatialIndex k{0}; k < getNz(); ++k) {
    auto ez_z_offset{(k + 0.5) * dz};
    auto hx_z_offset{(k + 0.5) * dz};
    for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
      auto ez_x_offset{i * dx};
      auto hx_x_offset{i * dx};
      _k_dot_r0_ez_yp(k, i) =
          (_k(0) * (max_x + ez_x_offset) + _k(1) * (max_y + ez_y_offset) +
           _k(2) * (max_z + ez_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hx_yp(k, i) =
          (_k(0) * (max_x + hx_x_offset) + _k(1) * (max_y + hx_y_offset) +
           _k(2) * (max_z + hx_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
    auto ex_z_offset{k * dz};
    auto hz_z_offset{k * dz};
    for (SpatialIndex i{0}; i < getNx(); ++i) {
      auto ex_x_offset{(i + 0.5) * dx};
      auto hz_x_offset{(i + 0.5) * dx};
      _k_dot_r0_ex_yp(k, i) = _k(0) * (max_x + ex_x_offset) +
                              _k(1) * (max_y + ex_y_offset) +
                              _k(2) * (max_z + ex_z_offset) - _l_0;
      _k_dot_r0_hz_yp(k, i) = _k(0) * (max_x + hz_x_offset) +
                              _k(1) * (max_y + hz_y_offset) +
                              _k(2) * (max_z + hz_z_offset) - _l_0;
    }
  }
}

void ArhTFSF1D::caculateKDotRZN(double dx, double dy, double dz) {
  if (getNx() == 1 && getNy() == 1) {
    // 1D
    auto min_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};
    _k_dot_r0_ex_zn(0, 0) = _k(2) * (min_z) / constant::C_0 - _l_0;
    _k_dot_r0_hy_zn(0, 0) = _k(2) * (min_z - 0.5 * dz) / constant::C_0 - _l_0;
    return;
  }

  auto min_x{_fdtd_domain_matrix(0, 0) + dx * getStartIndexX()};
  auto min_y{_fdtd_domain_matrix(0, 1) + dy * getStartIndexY()};
  auto min_z{_fdtd_domain_matrix(0, 2) + dz * getStartIndexZ()};
  double ex_z_offset{0.0};
  double ey_z_offset{0.0};
  double hx_z_offset{-0.5 * dz};
  double hy_z_offset{-0.5 * dz};

  for (SpatialIndex i{0}; i < getNx(); ++i) {
    auto ex_y_offset{(i + 0.5) * dy};
    auto hy_y_offset{(i + 0.5) * dy};
    for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
      auto ex_x_offset{j * dx};
      auto hy_x_offset{j * dx};
      _k_dot_r0_ex_zn(i, j) =
          (_k(0) * (min_x + ex_x_offset) + _k(1) * (min_y + ex_y_offset) +
           _k(2) * (min_z + ex_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hy_zn(i, j) =
          (_k(0) * (min_x + hy_x_offset) + _k(1) * (min_y + hy_y_offset) +
           _k(2) * (min_z + hy_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    auto ey_y_offset{i * dy};
    auto hx_y_offset{i * dy};
    for (SpatialIndex j{0}; j < getNy(); ++j) {
      auto ey_x_offset{(j + 0.5) * dx};
      auto hx_x_offset{(j + 0.5) * dx};
      _k_dot_r0_ey_zn(i, j) =
          (_k(0) * (min_x + ey_x_offset) + _k(1) * (min_y + ey_y_offset) +
           _k(2) * (min_z + ey_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hx_zn(i, j) =
          (_k(0) * (min_x + hx_x_offset) + _k(1) * (min_y + hx_y_offset) +
           _k(2) * (min_z + hx_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }
}

void ArhTFSF1D::caculateKDotRZP(double dx, double dy, double dz) {
  auto max_x{_fdtd_domain_matrix(0, 0) + dx * getStartIndexX()};
  auto max_y{_fdtd_domain_matrix(0, 1) + dy * getStartIndexY()};
  auto max_z{_fdtd_domain_matrix(0, 2) + dz * getEndIndexZ()};

  if (getNx() == 1 && getNy() == 1) {
    // 1D
    _k_dot_r0_ex_zp(0, 0) = _k(2) * (max_z) / constant::C_0 - _l_0;
    _k_dot_r0_hy_zp(0, 0) = _k(2) * (max_z + 0.5 * dz) / constant::C_0 - _l_0;
  }

  double ex_z_offset{0.0};
  double ey_z_offset{0.0};
  double hx_z_offset{0.5 * dz};
  double hy_z_offset{0.5 * dz};

  for (SpatialIndex i{0}; i < getNx(); ++i) {
    auto ex_y_offset{(i + 0.5) * dy};
    auto hy_y_offset{(i + 0.5) * dy};
    for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
      auto ex_x_offset{j * dx};
      auto hy_x_offset{j * dx};
      _k_dot_r0_ex_zp(i, j) =
          (_k(0) * (max_x + ex_x_offset) + _k(1) * (max_y + ex_y_offset) +
           _k(2) * (max_z + ex_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hy_zp(i, j) =
          (_k(0) * (max_x + hy_x_offset) + _k(1) * (max_y + hy_y_offset) +
           _k(2) * (max_z + hy_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }

  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    auto ey_y_offset{i * dy};
    auto hx_y_offset{i * dy};
    for (SpatialIndex j{0}; j < getNy(); ++j) {
      auto ey_x_offset{(j + 0.5) * dx};
      auto hx_x_offset{(j + 0.5) * dx};
      _k_dot_r0_ey_zp(i, j) =
          (_k(0) * (max_x + ey_x_offset) + _k(1) * (max_y + ey_y_offset) +
           _k(2) * (max_z + ey_z_offset)) /
              constant::C_0 -
          _l_0;
      _k_dot_r0_hx_zp(i, j) =
          (_k(0) * (max_x + hx_x_offset) + _k(1) * (max_y + hx_y_offset) +
           _k(2) * (max_z + hx_z_offset)) /
              constant::C_0 -
          _l_0;
    }
  }
}

void ArhTFSF1D::updateEiX(double time) {
  for (SpatialIndex j{0}; j < getNy(); ++j) {
    for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
      _eyi_xn(j, k) =
          _ey_i0 * _waveform->getValueByTime(time - _k_dot_r0_ey_xn(j, k));
      _eyi_xp(j, k) =
          _ey_i0 * _waveform->getValueByTime(time - _k_dot_r0_ey_xp(j, k));
    }
  }
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    for (SpatialIndex k{0}; k < getNz(); ++k) {
      _ezi_xn(j, k) =
          _ez_i0 * _waveform->getValueByTime(time - _k_dot_r0_ez_xn(j, k));
      _ezi_xp(j, k) =
          _ez_i0 * _waveform->getValueByTime(time - _k_dot_r0_ez_xp(j, k));
    }
  }
}

void ArhTFSF1D::updateEiY(double time) {
  for (SpatialIndex k{0}; k < getNz(); ++k) {
    for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
      _ezi_yn(k, i) =
          _ez_i0 * _waveform->getValueByTime(time - _k_dot_r0_ez_yn(k, i));
      _ezi_yp(k, i) =
          _ez_i0 * _waveform->getValueByTime(time - _k_dot_r0_ez_yp(k, i));
    }
  }

  for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
    for (SpatialIndex i{0}; i < getNx(); ++i) {
      _exi_yn(k, i) =
          _ex_i0 * _waveform->getValueByTime(time - _k_dot_r0_ex_yn(k, i));
      _exi_yp(k, i) =
          _ex_i0 * _waveform->getValueByTime(time - _k_dot_r0_ex_yp(k, i));
    }
  }
}

void ArhTFSF1D::updateEiZ(double time) {
  for (SpatialIndex i{0}; i < getNx(); ++i) {
    for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
      _exi_zn(i, j) =
          _ex_i0 * _waveform->getValueByTime(time - _k_dot_r0_ex_zn(i, j));
      _exi_zp(i, j) =
          _ex_i0 * _waveform->getValueByTime(time - _k_dot_r0_ex_zp(i, j));
    }
  }
  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    for (SpatialIndex j{0}; j < getNy(); ++j) {
      _eyi_zn(i, j) =
          _ey_i0 * _waveform->getValueByTime(time - _k_dot_r0_ey_zn(i, j));
      _eyi_zp(i, j) =
          _ey_i0 * _waveform->getValueByTime(time - _k_dot_r0_ey_zp(i, j));
    }
  }
}

void ArhTFSF1D::updateHiX(double time) {
  for (SpatialIndex j{0}; j < getNy(); ++j) {
    for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
      _hzi_xn(j, k) =
          _hz_i0 * _waveform->getValueByTime(time - _k_dot_r0_hz_xn(j, k));
      _hzi_xp(j, k) =
          _hz_i0 * _waveform->getValueByTime(time - _k_dot_r0_hz_xp(j, k));
    }
  }
  for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
    for (SpatialIndex k{0}; k < getNz(); ++k) {
      _hyi_xn(j, k) =
          _hy_i0 * _waveform->getValueByTime(time - _k_dot_r0_hy_xn(j, k));
      _hyi_xp(j, k) =
          _hy_i0 * _waveform->getValueByTime(time - _k_dot_r0_hy_xp(j, k));
    }
  }
}

void ArhTFSF1D::updateHiY(double time) {
  for (SpatialIndex k{0}; k < getNz(); ++k) {
    for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
      _hxi_yn(k, i) =
          _hx_i0 * _waveform->getValueByTime(time - _k_dot_r0_hx_yn(k, i));
      _hxi_yp(k, i) =
          _hx_i0 * _waveform->getValueByTime(time - _k_dot_r0_hx_yp(k, i));
    }
  }

  for (SpatialIndex k{0}; k < getNz() + 1; ++k) {
    for (SpatialIndex i{0}; i < getNx(); ++i) {
      _hzi_yn(k, i) =
          _hz_i0 * _waveform->getValueByTime(time - _k_dot_r0_hz_yn(k, i));
      _hzi_yp(k, i) =
          _hz_i0 * _waveform->getValueByTime(time - _k_dot_r0_hz_yp(k, i));
    }
  }
}

void ArhTFSF1D::updateHiZ(double time) {
  for (SpatialIndex i{0}; i < getNx(); ++i) {
    for (SpatialIndex j{0}; j < getNy() + 1; ++j) {
      _hyi_zn(i, j) =
          _hy_i0 * _waveform->getValueByTime(time - _k_dot_r0_hy_zn(i, j));
      _hyi_zp(i, j) =
          _hy_i0 * _waveform->getValueByTime(time - _k_dot_r0_hy_zp(i, j));
    }
  }
  for (SpatialIndex i{0}; i < getNx() + 1; ++i) {
    for (SpatialIndex j{0}; j < getNy(); ++j) {
      _hxi_zn(i, j) =
          _hx_i0 * _waveform->getValueByTime(time - _k_dot_r0_hx_zn(i, j));
      _hxi_zp(i, j) =
          _hx_i0 * _waveform->getValueByTime(time - _k_dot_r0_hx_zp(i, j));
    }
  }
}

}  // namespace xfdtd
