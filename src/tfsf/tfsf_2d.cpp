#include "tfsf/tfsf_2d.h"

#include <cmath>

#include "electromagnetic.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {
TFSF2D::TFSF2D(SpatialIndex distance_x, SpatialIndex distance_y, double phi_inc,
               double ez_i0, std::unique_ptr<Waveform> waveform)
    : TFSF(distance_x, distance_y, 0, constant::PI / 2, phi_inc, ez_i0, 0,
           std::move(waveform)) {}

void TFSF2D::init(const Cube *simulation_box, double dx, double dy, double dz,
                  double dt, TFSFBoundaryIndex tfsf_boundary_index) {
  defaultInitTFSF(simulation_box, dx, dy, dz, dt, tfsf_boundary_index);
  // TODO(franzero)
  _ceihi = -(getDt() / (constant::EPSILON_0 * getDx()));
  _chiei = -(getDt() / (constant::MU_0 * getDx()));

  // IMPORTANT: assume that dx is equal to dy.
  auto incident_phi{getIncidentPhi()};
  if (isGreaterOrEqual(incident_phi, 0.0, constant::TOLERABLE_EPSILON) &&
      isLessOrEqual(incident_phi, constant::PI / 2,
                    constant::TOLERABLE_EPSILON)) {
    _strike_index_x = getStartIndexX();
    _strike_index_y = getStartIndexY();
  } else if (isGreaterOrEqual(incident_phi, constant::PI / 2,
                              constant::TOLERABLE_EPSILON) &&
             isLessOrEqual(incident_phi, constant::PI,
                           constant::TOLERABLE_EPSILON)) {
    _strike_index_x = getEndIndexX();
    _strike_index_y = getStartIndexY();
  } else if (isGreaterOrEqual(incident_phi, constant::PI,
                              constant::TOLERABLE_EPSILON) &&
             isLessOrEqual(incident_phi, constant::PI * 1.5,
                           constant::TOLERABLE_EPSILON)) {
    _strike_index_x = getEndIndexX();
    _strike_index_y = getEndIndexY();
  } else {
    _strike_index_x = getStartIndexX();
    _strike_index_y = getEndIndexY();
  }

  auto diagonal_length{std::sqrt(std::pow(getNx(), 2) + std::pow(getNy(), 2))};
  // why is the length times 2 plus 4 and plus 1?
  // The number of E points on the diagonal = number of points on on side of %
  // TF / SF boundary. Add 4 to this, 2 on either end
  // TODO(franzero): Can observe reflection obviously while using Mur abosrbing
  // Bbundary bondition. A temporary solution is to extend the length.
  size_t temp_times{4};
  _auxiliary_array_size = std::ceil(diagonal_length * temp_times) + 4 + 1;
  _e_inc.resize(_auxiliary_array_size);
  _h_inc.resize(_auxiliary_array_size - 1);

  _projection_index_ez_xn.resize(getNy() + 1);
  _projection_index_hy_xn.resize(getNy() + 1);

  _projection_index_ez_yn.resize(getNx() + 1);
  _projection_index_hx_yn.resize(getNx() + 1);

  _projection_index_ez_xp.resize(getNy() + 1);
  _projection_index_hy_xp.resize(getNy() + 1);

  _projection_index_ez_yp.resize(getNx() + 1);
  _projection_index_hx_yp.resize(getNx() + 1);

  _weight_ez_xn.resize(getNy() + 1);
  _weight_hy_xn.resize(getNy() + 1);

  _weight_ez_yn.resize(getNx() + 1);
  _weight_hx_yn.resize(getNx() + 1);

  _weight_ez_xp.resize(getNy() + 1);
  _weight_hy_xp.resize(getNy() + 1);

  _weight_ez_yp.resize(getNx() + 1);
  _weight_hx_yp.resize(getNx() + 1);

  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};
  auto k_inc{getKVector()};

  for (auto i{li}; i < ri + 1; ++i) {
    auto vec_ez_yn{
        PointVector{i - getStrikeIndexX(), lj - getStrikeIndexY(), 0}};
    auto k_dot_ez_yn{k_inc.transpose() * vec_ez_yn};
    _projection_index_ez_yn[i - li] =
        std::floor(k_dot_ez_yn) + 2;  // 2 for boundary
    _weight_ez_yn[i - li] = k_dot_ez_yn - _projection_index_ez_yn[i - li] + 2;

    auto vec_hx_yn{
        PointVector{i - getStrikeIndexX(), lj - 0.5 - getStrikeIndexY(), 0}};
    auto k_dot_hx_yn{k_inc.transpose() * vec_hx_yn};
    _projection_index_hx_yn[i - li] = std::floor(k_dot_hx_yn) + 2;
    _weight_hx_yn[i - li] = k_dot_hx_yn - _projection_index_hx_yn[i - li] + 2;

    auto vec_ez_yp{
        PointVector{i - getStrikeIndexX(), rj - getStrikeIndexY(), 0}};
    auto k_dot_ez_yp{k_inc.transpose() * vec_ez_yp};
    _projection_index_ez_yp[i - li] = std::floor(k_dot_ez_yp) + 2;
    _weight_ez_yp[i - li] = k_dot_ez_yp - _projection_index_ez_yp[i - li] + 2;

    auto vec_hx_yp{
        PointVector{i - getStrikeIndexX(), rj + 0.5 - getStrikeIndexY(), 0}};
    auto k_dot_hx_yp{k_inc.transpose() * vec_hx_yp};
    _projection_index_hx_yp[i - li] = std::floor(k_dot_hx_yp) + 2;
    _weight_hx_yp[i - li] = k_dot_hx_yp - _projection_index_hx_yp[i - li] + 2;
  }

  for (auto j{lj}; j < rj + 1; ++j) {
    auto vec_ez_xn{
        PointVector{li - getStrikeIndexX(), j - getStrikeIndexY(), 0}};
    auto k_dot_ez_xn{k_inc.transpose() * vec_ez_xn};
    _projection_index_ez_xn[j - lj] = std::floor(k_dot_ez_xn) + 2;
    _weight_ez_xn[j - lj] = k_dot_ez_xn - _projection_index_ez_xn[j - lj] + 2;

    auto vec_hy_xn{
        PointVector{li - 0.5 - getStrikeIndexX(), j - getStrikeIndexY(), 0}};
    auto k_dot_hy_xn{k_inc.transpose() * vec_hy_xn};
    _projection_index_hy_xn[j - lj] = std::floor(k_dot_hy_xn) + 2;
    _weight_hy_xn[j - lj] = k_dot_hy_xn - _projection_index_hy_xn[j - lj] + 2;

    auto vec_ez_xp{
        PointVector{ri - getStrikeIndexX(), j - getStrikeIndexY(), 0}};
    auto k_dot_ez_xp{k_inc.transpose() * vec_ez_xp};
    _projection_index_ez_xp[j - lj] = std::floor(k_dot_ez_xp) + 2;
    _weight_ez_xp[j - lj] = k_dot_ez_xp - _projection_index_ez_xp[j - lj] + 2;

    auto vec_hy_xp{
        PointVector{ri + 0.5 - getStrikeIndexX(), j - getStrikeIndexY(), 0}};
    auto k_dot_hy_xp{k_inc.transpose() * vec_hy_xp};
    _projection_index_hy_xp[j - lj] = std::floor(k_dot_hy_xp) + 2;
    _weight_hy_xp[j - lj] = k_dot_hy_xp - _projection_index_hy_xp[j - lj] + 2;
  }
}

void TFSF2D::updateIncidentField(size_t current_time_step) {
  _e_inc[0] = getIncidentFieldWaveformValueByTime(current_time_step * getDt());
  // 1D Mur Absorbing Boundary Condition
  auto alpha{constant::C_0 * getDt() - getDx() / constant::C_0 * getDt() +
             getDx()};
  auto temp{_e_inc[_e_inc.size() - 2]};
  for (auto i{1}; i < _e_inc.size() - 1; ++i) {
    _e_inc[i] = _ceie * _e_inc[i] + _ceihi * (_h_inc[i] - _h_inc[i - 1]);
  }
  _e_inc[_e_inc.size() - 1] =
      temp + alpha * (_e_inc[_e_inc.size() - 2] - _e_inc[_e_inc.size() - 1]);
  for (auto i{0}; i < _h_inc.size(); ++i) {
    _h_inc[i] = _chih * _h_inc[i] + _chiei * (_e_inc[i + 1] - _e_inc[i]);
  }
}

void TFSF2D::updateH() {
  auto li{getStartIndexX()};
  auto ri{getEndIndexX()};
  auto lj{getStartIndexY()};
  auto rj{getEndIndexY()};

  for (auto i{li}; i < ri + 1; ++i) {
    auto p{_projection_index_ez_yn[i - li]};
    auto w{_weight_ez_yn[i - li]};
    auto ezi_yn{(1 - w) * _e_inc[p] + w * _e_inc[p + 1]};
    auto &hx_yn{getHx(i, lj - 1, 0)};
    hx_yn += (getDt() / (constant::MU_0 * getDy())) * ezi_yn;

    p = _projection_index_ez_yp[i - li];
    w = _weight_ez_yp[i - li];
    auto ezi_yp{(1 - w) * _e_inc[p] + w * _e_inc[p + 1]};
    auto &hx_yp{getHx(i, rj, 0)};
    hx_yp -= (getDt() / (constant::MU_0 * getDy())) * ezi_yp;
  }

  for (auto j{lj}; j < rj + 1; ++j) {
    auto p{_projection_index_hy_xn[j - lj]};
    auto w{_weight_hy_xn[j - lj]};
    auto ezi_xn{(1 - w) * _e_inc[p] + w * _e_inc[p + 1]};
    auto &hy_xn{getHy(li - 1, j, 0)};
    hy_xn -= (getDt() / (constant::MU_0 * getDx())) * ezi_xn;

    p = _projection_index_hy_xp[j - lj];
    w = _weight_hy_xp[j - lj];
    auto ezi_xp{(1 - w) * _e_inc[p] + w * _e_inc[p + 1]};
    auto &hy_xp{getHy(ri, j, 0)};
    hy_xp += (getDt() / (constant::MU_0 * getDx())) * ezi_xp;
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
  auto sin_phi{std::sin(getIncidentPhi())};
  auto cos_phi{std::cos(getIncidentPhi())};

  for (auto i{li}; i < ri + 1; ++i) {
    auto p{_projection_index_ez_yn[i - li]};
    auto w{_weight_ez_yn[i - li]};
    auto hxi_yn{(1 - w) * sin_phi * _h_inc[p] + w * sin_phi * _h_inc[p + 1]};
    auto &ex_yn{getEz(i, lj, 0)};
    ex_yn += (dt / (constant::EPSILON_0 * dy)) * hxi_yn;

    p = _projection_index_ez_yp[i - li];
    w = _weight_ez_yp[i - li];
    auto hxi_yp{(1 - w) * sin_phi * _h_inc[p] + w * sin_phi * _h_inc[p + 1]};
    auto &ex_yp{getEz(i, rj, 0)};
    ex_yp -= (dt / (constant::EPSILON_0 * dy)) * hxi_yp;
  }

  for (auto j{lj}; j < rj + 1; ++j) {
    auto p{_projection_index_ez_xn[j - lj]};
    auto w{_weight_ez_xn[j - lj]};
    auto hyi_xn{-(1 - w) * cos_phi * _h_inc[p] - w * cos_phi * _h_inc[p + 1]};
    auto &ez_xn{getEz(li, j, 0)};
    ez_xn -= (dt / (constant::EPSILON_0 * dx)) * hyi_xn;

    p = _projection_index_ez_xp[j - lj];
    w = _weight_ez_xp[j - lj];
    auto hyi_xp{-(1 - w) * cos_phi * _h_inc[p] - w * cos_phi * _h_inc[p + 1]};
    auto &ez_xp{getEz(ri, j, 0)};
    ez_xp += (dt / (constant::EPSILON_0 * dx)) * hyi_xp;
  }
}
}  // namespace xfdtd