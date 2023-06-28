#include "tfsf/tfsf_2d.h"

#include <fstream>
#include <iostream>
#include <istream>

#include "util/constant.h"
namespace xfdtd {
TFSF2D::TFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
               double ez_i0, std::unique_ptr<Waveform> waveform)
    : TFSF(distance_x, distance_y, 0, ez_i0, constant::PI / 2, phi_inc, 0,
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
  auto k_e{PointVector{sin_psi, cos_psi}};
  _transform_e = xt::linalg::dot(u, k_e);
  _transform_h = xt::linalg::cross(getKVector(), _transform_e);
}

void TFSF2D::init(double dx, double dy, double dz, double dt,
                  std::unique_ptr<GridBox> tfsf_grid_box) {
  defaultInitTFSF(dx, dy, dz, dt, std::move(tfsf_grid_box));
  const auto nx{getNx()};
  const auto ny{getNy()};

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
      static_cast<size_t>(std::ceil(diagonal_length)) + 4 + 1;
  _e_inc.resize({_auxiliary_array_size});
  _h_inc.resize({_auxiliary_array_size - 1});

  _ez_inc.resize({_auxiliary_array_size});
  _hx_inc.resize({_auxiliary_array_size - 1});
  _hy_inc.resize({_auxiliary_array_size - 1});
  auto extra_distance{2 * getKInc() / ratio_delta};
  PointVector origination_point;
  if (k_inc(0) >= 0 && k_inc(1) >= 0) {
    // xp yp
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getStartIndexY()), 0};

  } else if (k_inc(0) >= 0 && k_inc(1) < 0) {
    // xp yn
    origination_point = PointVector{static_cast<double>(getStartIndexX()),
                                    static_cast<double>(getEndIndexY()), 0};
  } else if (k_inc(0) < 0 && k_inc(1) >= 0) {
    // xn yp
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getStartIndexY()), 0};
  } else if (k_inc(0) < 0 && k_inc(1) < 0) {
    // xn yn
    origination_point = PointVector{static_cast<double>(getEndIndexX()),
                                    static_cast<double>(getEndIndexY()), 0};
  }
  origination_point -= extra_distance;

  _projection_x_full = xt::zeros<double>({nx + 1});
  _projection_y_full = xt::zeros<double>({ny + 1});
  _projection_x_half = xt::zeros<double>({nx + 2});
  _projection_y_half = xt::zeros<double>({ny + 2});
  const auto k_inc_x{k_inc[0]};
  const auto k_inc_y{k_inc[1]};

  for (size_t i{0}; i < nx + 1; ++i) {
    _projection_x_full(i) =
        (i + getStartIndexX() - origination_point(0)) * k_inc_x * ratio_delta;
  }
  for (size_t i{0}; i < ny + 1; ++i) {
    _projection_y_full(i) =
        (i + getStartIndexY() - origination_point(1)) * k_inc_y * ratio_delta;
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
  // ez
  for (auto i{getStartIndexX()}; i < getEndIndexX(); ++i) {
    getIncidentEz(i, getStartIndexY(), 0);
  }
  for (auto i{getStartIndexX()}; i < getEndIndexX(); ++i) {
    getIncidentEz(i, getEndIndexY(), 0);
  }

  for (SpatialIndex i{getStartIndexY()}; i < getEndIndexY() + 1; ++i) {
    getIncidentHy(getStartIndexX() - 1, i, 0);
  }
  std::cout << std::endl;

  _scaled_dl = getDx() / ratio_delta;
  _ceihi = -(dt / (constant::EPSILON_0 * _scaled_dl));
  _chiei = -(dt / (constant::MU_0 * _scaled_dl));
}

void TFSF2D::updateIncidentField(size_t current_time_step) {
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

  _ez_inc = _transform_e(2) * _e_inc;
  _hx_inc = _transform_h(0) * _h_inc;
  _hy_inc = _transform_h(1) * _h_inc;
}

void TFSF2D::updateH() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};

  auto dt{getDt()};
  auto dx{getDx()};
  auto dy{getDy()};

  auto cbx{dt / (constant::MU_0 * dx)};
  auto cby{dt / (constant::MU_0 * dy)};

  // xn
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    auto& hy_xn{getHy(li - 1, j, 0)};
    auto ezi{getIncidentEz(li, j, 0)};
    hy_xn -= cbx * ezi;
  }

  // xp
  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    auto& hy_xp{getHy(ri, j, 0)};
    auto ezi{getIncidentEz(ri, j, 0)};
    hy_xp += cbx * ezi;
  }

  // yn
  for (auto i{li}; i < ri + 1; ++i) {
    auto& hx_yn{getHx(i, lj - 1, 0)};
    auto ezi{getIncidentEz(i, lj, 0)};
    hx_yn += cby * ezi;
  }

  // yp
  for (auto i{li}; i < ri + 1; ++i) {
    auto& hx_yp{getHx(i, rj, 0)};
    auto ezi{getIncidentEz(i, rj, 0)};
    hx_yp -= cby * ezi;
  }
}

void TFSF2D::updateE() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};

  auto dt{getDt()};
  auto dx{getDx()};
  auto dy{getDy()};

  auto cax{dt / (constant::EPSILON_0 * dx)};
  auto cay{dt / (constant::EPSILON_0 * dy)};

  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    auto& ez_xn{getEz(li, j, 0)};
    auto hyi{getIncidentHy(li - 1, j, 0)};
    ez_xn -= cax * hyi;
  }

  for (SpatialIndex j{lj}; j < rj + 1; ++j) {
    auto& ez_xp{getEz(ri, j, 0)};
    auto hyi{getIncidentHy(ri, j, 0)};
    ez_xp += cax * hyi;
  }

  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    auto& ez_yn{getEz(i, lj, 0)};
    auto hxi{getIncidentHx(i, lj - 1, 0)};
    ez_yn += cay * hxi;
  }

  for (SpatialIndex i{li}; i < ri + 1; ++i) {
    auto& ez_yp{getEz(i, rj, 0)};
    auto hxi{getIncidentHx(i, rj, 0)};
    ez_yp -= cay * hxi;
  }
}

double TFSF2D::getIncidentEz(int i, int j, int k) {
  // 0 0 0.5
  i = i - getStartIndexX();
  j = j - getStartIndexY();
  auto projection{_projection_x_full(i) + _projection_y_full(j)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _ez_inc(index) + weight * _ez_inc(index + 1);
}

double TFSF2D::getIncidentHx(int i, int j, int k) {
  // 0 0.5 0.5
  i = i - getStartIndexX();
  j = j - getStartIndexY() + 1;
  auto projection{_projection_x_full(i) + _projection_y_half(j)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hx_inc(index) + weight * _hx_inc(index + 1);
}

double TFSF2D::getIncidentHy(int i, int j, int k) {
  // 0.5 0 0.5
  i = i - getStartIndexX() + 1;
  j = j - getStartIndexY();
  auto projection{_projection_x_half(i) + _projection_y_full(j)};
  auto index{static_cast<SpatialIndex>(projection)};
  auto weight{projection - index};
  return (1 - weight) * _hy_inc(index) + weight * _hy_inc(index + 1);
}

}  // namespace xfdtd