#include "tfsf/tfsf_3d.h"

#include <cmath>
#include <utility>
#include <xtensor/xfixed.hpp>

#include "tfsf/tfsf.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"
#include "xtensor-blas/xlinalg.hpp"

namespace xfdtd {

TFSF3D::TFSF3D(SpatialIndex distance_x, SpatialIndex distance_y,
               SpatialIndex distance_z, double e_0, double theta_inc,
               double phi_inc, double psi, std::unique_ptr<Waveform> waveform)
    : TFSF(distance_x, distance_y, distance_z, e_0, theta_inc, phi_inc, psi,
           std::move(waveform)) {
  auto sin_phi{getIncidentSinPhi()};
  auto cos_phi{getIncidentCosPhi()};
  auto sin_theta{getIncidentSinTheta()};
  auto cos_theta{getIncidentCosTheta()};
  xt::xtensor<double, 2> u{{-sin_phi, cos_theta * cos_phi, sin_theta * cos_phi},
                           {cos_phi, cos_theta * sin_phi, sin_theta * sin_phi},
                           {0, -sin_theta, cos_theta}};

  auto sin_psi{getIncidentSinPsi()};
  auto cos_psi{getIncidentCosPsi()};
  auto k_e{PointVector{sin_psi, cos_psi, 0}};
  _transform_e = xt::linalg::dot(u, k_e);
  _transform_h = xt::linalg::cross(getKVector(), _transform_e);
};

void TFSF3D::init(double dx, double dy, double dz, double dt,
                  std::unique_ptr<GridBox> tfsf_grid_box) {
  defaultInitTFSF(dx, dy, dz, dt, std::move(tfsf_grid_box));
  const auto nx{getNx()};
  const auto ny{getNy()};
  const auto nz{getNz()};

  // IMPORTANT: assume that dx is equal to dy.
  const auto incident_phi{getIncidentPhi()};
  const auto incident_theta{getIncidentTheta()};
  const auto ratio_delta{
      1 / sqrt(pow(sin(incident_theta), 4) *
                   (pow(cos(incident_phi), 4) + pow(sin(incident_phi), 4)) +
               pow(cos(incident_theta), 4))};
  const auto k_inc{getKInc()};

  const auto diagonal_length{
      sqrt(pow(getNx(), 2) + pow(getNy(), 2) + sqrt(pow(getNz(), 2))) *
      ratio_delta};
  _auxiliary_array_size =
      static_cast<size_t>(std::ceil(diagonal_length)) * 4 + 4 + 1;
  _e_inc.resize({_auxiliary_array_size});
  _ex_inc.resize({_auxiliary_array_size});
  _ey_inc.resize({_auxiliary_array_size});
  _ez_inc.resize({_auxiliary_array_size});
  _h_inc.resize({_auxiliary_array_size - 1});
  _hx_inc.resize({_auxiliary_array_size - 1});
  _hy_inc.resize({_auxiliary_array_size - 1});
  _hz_inc.resize({_auxiliary_array_size - 1});
  auto extra_distance{2 * getKInc() / ratio_delta};
  PointVector origination_point;
  if (k_inc(0) >= 0 && k_inc(1) >= 0 && k_inc(2) >= 0) {
    // xp yp zp
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getStartIndexY()),
                                    static_cast<double>(getStartIndexZ())};
  } else if (k_inc(0) >= 0 && k_inc(1) >= 0 && k_inc(2) < 0) {
    // xp yp zn
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getStartIndexY()),
                                    static_cast<double>(getEndIndexZ())};
  } else if (k_inc(0) >= 0 && k_inc(1) < 0 && k_inc(2) >= 0) {
    // xp yn zp
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getEndIndexY()),
                                    static_cast<double>(getStartIndexZ())};
  } else if (k_inc(0) >= 0 && k_inc(1) < 0 && k_inc(2) < 0) {
    // xp yn zn
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getEndIndexY()),
                                    static_cast<double>(getEndIndexZ())};
  } else if (k_inc(0) < 0 && k_inc(1) >= 0 && k_inc(2) >= 0) {
    // xn yp zp
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getStartIndexY()),
                                    static_cast<double>(getStartIndexZ())};
  } else if (k_inc(0) < 0 && k_inc(1) >= 0 && k_inc(2) < 0) {
    // xn yp zn
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getStartIndexY()),
                                    static_cast<double>(getEndIndexZ())};
  } else if (k_inc(0) < 0 && k_inc(1) < 0 && k_inc(2) >= 0) {
    // xn yn zp
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getEndIndexY()),
                                    static_cast<double>(getStartIndexZ())};
  } else {
    // xn yn zn
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getEndIndexY()),
                                    static_cast<double>(getEndIndexZ())};
  }
  origination_point -= extra_distance;

  _projection_x_full = xt::zeros<double>({nx + 1});
  _projection_y_full = xt::zeros<double>({ny + 1});
  _projection_z_full = xt::zeros<double>({nz + 1});
  _projection_x_half = xt::zeros<double>({nx + 2});
  _projection_y_half = xt::zeros<double>({ny + 2});
  _projection_z_half = xt::zeros<double>({nz + 2});
  const auto k_inc_x{k_inc[0]};
  const auto k_inc_y{k_inc[1]};
  const auto k_inc_z{k_inc[2]};
  for (size_t i{0}; i < nx + 1; ++i) {
    _projection_x_full(i) =
        (i + getStartIndexX() - origination_point(0)) * k_inc_x * ratio_delta;
  }
  for (size_t i{0}; i < ny + 1; ++i) {
    _projection_y_full(i) =
        (i + getStartIndexY() - origination_point(1)) * k_inc_y * ratio_delta;
  }
  for (size_t i{0}; i < nz + 1; ++i) {
    _projection_z_full(i) =
        (i + getStartIndexZ() - origination_point(2)) * k_inc_z * ratio_delta;
  }
  for (size_t i{0}; i < nx + 2; ++i) {
    _projection_x_half(i) =
        (i + getStartIndexX() - origination_point(0) - 0.5) * k_inc_x *
        ratio_delta;
  }
  for (size_t i{0}; i < ny + 2; ++i) {
    _projection_y_half(i) =
        (i + getStartIndexY() - origination_point(1) - 0.5) * k_inc_y *
        ratio_delta;
  }
  for (size_t i{0}; i < nz + 2; ++i) {
    _projection_z_half(i) =
        (i + getStartIndexZ() - origination_point(2) - 0.5) * k_inc_z *
        ratio_delta;
  }

  _scaled_dl = getDx() / ratio_delta;
  _ceihi = -(dt / (constant::EPSILON_0 * _scaled_dl));
  _chiei = -(dt / (constant::MU_0 * _scaled_dl));
}

void TFSF3D::updateIncidentField(size_t current_time_step) {
  auto dt{getDt()};
  _e_inc[0] = getIncidentFieldWaveformValueByTime(current_time_step * getDt());
  // 1D Mur Absorbing Boundary Condition
  for (auto i{1}; i < _e_inc.size() - 1; ++i) {
    _e_inc[i] = _ceie * _e_inc[i] + _ceihi * (_h_inc[i] - _h_inc[i - 1]);
  }
  _e_inc[_e_inc.size() - 1] =
      _e_inc[_e_inc.size() - 1] -
      (constant::C_0 * dt / _scaled_dl) *
          (_e_inc[_e_inc.size() - 1] - _e_inc[_e_inc.size() - 2]);

  for (auto i{0}; i < _h_inc.size(); ++i) {
    _h_inc[i] = _chih * _h_inc[i] + _chiei * (_e_inc[i + 1] - _e_inc[i]);
  }

  _ex_inc = _transform_e(0) * _e_inc;
  _ey_inc = _transform_e(1) * _e_inc;
  _ez_inc = _transform_e(2) * _e_inc;
  _hx_inc = _transform_h(0) * _h_inc;
  _hy_inc = _transform_h(1) * _h_inc;
  _hz_inc = _transform_h(2) * _h_inc;
}

void TFSF3D::updateH() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};
  auto lk{getStartIndexZ()};
  auto rk{getEndIndexZ()};

  auto dt{getDt()};
  auto dx{getDx()};
  auto dy{getDy()};
  auto dz{getDz()};

  auto cbx{dt / (constant::MU_0 * dx)};
  auto cby{dt / (constant::MU_0 * dy)};
  auto cbz{dt / (constant::MU_0 * dz)};

  // xn
  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& hz_xn{getHz(li - 1, j, k)};
      auto eyi{getIncidentEy(li, j, k)};
      hz_xn += cbx * eyi;
    }
  }
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& hy_xn{getHy(li - 1, j, k)};
      auto ezi{getIncidentEz(li, j, k)};
      hy_xn -= cbx * ezi;
    }
  }

  // xp
  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& hz_xp{getHz(ri, j, k)};
      auto eyi{getIncidentEy(ri, j, k)};
      hz_xp -= cbx * eyi;
    }
  }
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& hy_xp{getHy(ri, j, k)};
      auto ezi{getIncidentEz(ri, j, k)};
      hy_xp += cbx * ezi;
    }
  }

  // yn
  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri + 1; ++i) {
      auto& hx_yn{getHx(i, lj - 1, k)};
      auto ezi{getIncidentEz(i, lj, k)};
      hx_yn += cby * ezi;
    }
  }
  for (auto k{lk}; k < rk + 1; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto& hz_yn{getHz(i, lj - 1, k)};
      auto exi{getIncidentEx(i, lj, k)};
      hz_yn -= cby * exi;
    }
  }

  // yp
  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri + 1; ++i) {
      auto& hx_yp{getHx(i, rj, k)};
      auto ezi{getIncidentEz(i, rj, k)};
      hx_yp -= cby * ezi;
    }
  }
  for (auto k{lk}; k < rk + 1; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto& hz_yp{getHz(i, rj, k)};
      auto exi{getIncidentEx(i, rj, k)};
      hz_yp += cby * exi;
    }
  }

  // zn
  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj + 1; ++j) {
      auto& hy_zn{getHy(i, j, lk - 1)};
      auto exi{getIncidentEx(i, j, lk)};
      hy_zn += cbz * exi;
    }
  }
  for (auto i{li}; i < ri + 1; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto& hx_zn{getHx(i, j, lk - 1)};
      auto eyi{getIncidentEy(i, j, lk)};
      hx_zn -= cbz * eyi;
    }
  }

  // zp
  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj + 1; ++j) {
      auto& hy_zp{getHy(i, j, rk)};
      auto exi{getIncidentEx(i, j, rk)};
      hy_zp -= cbz * exi;
    }
  }
  for (auto i{li}; i < ri + 1; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto& hx_zp{getHx(i, j, rk)};
      auto eyi{getIncidentEy(i, j, rk)};
      hx_zp += cbz * eyi;
    }
  }
}

void TFSF3D::updateE() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};
  auto lk{getStartIndexZ()};
  auto rk{getEndIndexZ()};

  auto dt{getDt()};
  auto dx{getDx()};
  auto dy{getDy()};
  auto dz{getDz()};

  auto cax{dt / (constant::EPSILON_0 * dx)};
  auto cay{dt / (constant::EPSILON_0 * dy)};
  auto caz{dt / (constant::EPSILON_0 * dz)};

  // xn
  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& ey_xn{getEy(li, j, k)};
      auto hzi{getIncidentHz(li - 1, j, k)};
      ey_xn += cax * hzi;
    }
  }
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& ez_xn{getEz(li, j, k)};
      auto hyi{getIncidentHy(li - 1, j, k)};
      ez_xn -= cax * hyi;
    }
  }

  // xp
  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (SpatialIndex k{lk}; k < rk + 1; ++k) {
      auto& ey_xp{getEy(ri, j, k)};
      auto hzi{getIncidentHz(ri, j, k)};
      ey_xp -= cax * hzi;
    }
  }
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    for (SpatialIndex k{lk}; k < rk; ++k) {
      auto& ez_xp{getEz(ri, j, k)};
      auto hyi{getIncidentHy(ri, j, k)};
      ez_xp += cax * hyi;
    }
  }

  // yn
  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri + 1; ++i) {
      auto& ez_yn{getEz(i, lj, k)};
      auto hxi{getIncidentHx(i, lj - 1, k)};
      ez_yn += cay * hxi;
    }
  }
  for (SpatialIndex k{lk}; k < rk + 1; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      auto& ex_yn{getEx(i, lj, k)};
      auto hzi{getIncidentHz(i, lj - 1, k)};
      ex_yn -= cay * hzi;
    }
  }

  // yp
  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri + 1; ++i) {
      auto& ez_yp{getEz(i, rj, k)};
      auto hxi{getIncidentHx(i, rj, k)};
      ez_yp -= cay * hxi;
    }
  }
  for (SpatialIndex k{lk}; k < rk + 1; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      auto& ex_yp{getEx(i, rj, k)};
      auto hzi{getIncidentHz(i, rj, k)};
      ex_yp += cay * hzi;
    }
  }

  // zn
  for (SpatialIndex i{li}; i < ri; ++i) {
    for (SpatialIndex j{lj}; j < rj + 1; ++j) {
      auto& ex_zn{getEx(i, j, lk)};
      auto hyi{getIncidentHy(i, j, lk - 1)};
      ex_zn += caz * hyi;
    }
  }
  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    for (SpatialIndex j{lj}; j < rj; ++j) {
      auto& ey_zn{getEy(i, j, lk)};
      auto hxi{getIncidentHx(i, j, lk - 1)};
      ey_zn -= caz * hxi;
    }
  }

  // zp
  for (SpatialIndex i{li}; i < ri; ++i) {
    for (SpatialIndex j{lj}; j < rj + 1; ++j) {
      auto& ex_zp{getEx(i, j, rk)};
      auto hyi{getIncidentHy(i, j, rk)};
      ex_zp -= caz * hyi;
    }
  }
  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    for (SpatialIndex j{lj}; j < rj; ++j) {
      auto& ey_zp{getEy(i, j, rk)};
      auto hxi{getIncidentHx(i, j, rk)};
      ey_zp += caz * hxi;
    }
  }
}

double TFSF3D::getIncidentEx(int i, int j, int k) {
  // 0.5 0 0
  i = i - getStartIndexX() + 1;
  j = j - getStartIndexY();
  k = k - getStartIndexZ();
  auto projection{_projection_x_half(i) + _projection_y_full(j) +
                  _projection_z_full(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  auto e_1{getIncidentCosPsi() * getIncidentCosTheta() * getIncidentCosPhi()};
  auto e_2{-getIncidentSinPsi() * getIncidentSinPhi()};
  auto ex{1};
  return (1 - weight) * _ex_inc(index) + weight * _ex_inc(index + 1);
}

double TFSF3D::getIncidentEy(int i, int j, int k) {
  // 0 0.5 0
  i = i - getStartIndexX();
  j = j - getStartIndexY() + 1;
  k = k - getStartIndexZ();
  auto projection{_projection_x_full(i) + _projection_y_half(j) +
                  _projection_z_full(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ey_inc(index) + weight * _ey_inc(index + 1);
}

double TFSF3D::getIncidentEz(int i, int j, int k) {
  // 0 0 0.5
  i = i - getStartIndexX();
  j = j - getStartIndexY();
  k = k - getStartIndexZ() + 1;
  auto projection{_projection_x_full(i) + _projection_y_full(j) +
                  _projection_z_half(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ez_inc(index) + weight * _ez_inc(index + 1);
}

double TFSF3D::getIncidentHx(int i, int j, int k) {
  // 0 0.5 0.5
  i = i - getStartIndexX();
  j = j - getStartIndexY() + 1;
  k = k - getStartIndexZ() + 1;
  auto projection{_projection_x_full(i) + _projection_y_half(j) +
                  _projection_z_half(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hx_inc(index) + weight * _hx_inc(index + 1);
}

double TFSF3D::getIncidentHy(int i, int j, int k) {
  // 0.5 0 0.5
  i = i - getStartIndexX() + 1;
  j = j - getStartIndexY();
  k = k - getStartIndexZ() + 1;
  auto projection{_projection_x_half(i) + _projection_y_full(j) +
                  _projection_z_half(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hy_inc(index) + weight * _hy_inc(index + 1);
}

double TFSF3D::getIncidentHz(int i, int j, int k) {
  // 0.5 0.5 0
  i = i - getStartIndexX() + 1;
  j = j - getStartIndexY() + 1;
  k = k - getStartIndexZ();
  auto projection{_projection_x_half(i) + _projection_y_half(j) +
                  _projection_z_full(k)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hz_inc(index) + weight * _hz_inc(index + 1);
}

}  // namespace xfdtd
