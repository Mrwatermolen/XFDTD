#include "nffft/nffft_fd.h"

#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <utility>
#include <xtensor/xadapt.hpp>
#include <xtensor/xnpy.hpp>

#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {

void NffftFd::init(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                   std::shared_ptr<const EMF> emf) {
  defaultInit(std::move(grid_space), std::move(fdtd_basic_coff),
              std::move(emf));

  _f_theta.resize({_number_frequencies, _number_theta, _number_phi});
  _f_phi.resize({_number_frequencies, _number_theta, _number_phi});
  _a_theta.resize({_number_frequencies, _number_theta, _number_phi});
  _a_phi.resize({_number_frequencies, _number_theta, _number_phi});
  _a_x.resize({_number_frequencies, _number_theta, _number_phi});
  _a_y.resize({_number_frequencies, _number_theta, _number_phi});
  _a_z.resize({_number_frequencies, _number_theta, _number_phi});
  _f_x.resize({_number_frequencies, _number_theta, _number_phi});
  _f_y.resize({_number_frequencies, _number_theta, _number_phi});
  _f_z.resize({_number_frequencies, _number_theta, _number_phi});

  _jy_xn.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jy_xp.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xn.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xp.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _jz_yn.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jz_yp.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jx_yn.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jx_yp.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});

  _jx_zn.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jx_zp.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jy_zn.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jy_zp.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});

  _my_xn.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _my_xp.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xn.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xp.resize({_number_frequencies, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _mz_yn.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mz_yp.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mx_yn.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mx_yp.resize({_number_frequencies, static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});

  _mx_zn.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _mx_zp.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _my_zn.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _my_zp.resize({_number_frequencies, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});

  initDFT();
}

void NffftFd::initDFT() {
  using namespace std::complex_literals;
  _frequency_transform_j.resize({_number_frequencies, getTotalTimeSteps()});
  _frequency_transform_m.resize({_number_frequencies, getTotalTimeSteps()});
  const auto dt{getDt()};
  for (int n{0}; n < _number_frequencies; ++n) {
    for (size_t t{0}; t < getTotalTimeSteps(); ++t) {
      // TODO(franzero): add or subtract 0.5?
      _frequency_transform_j(n, t) =
          dt * std::exp(-1i * 2.0 * constant::PI *
                        (_frequencies(n) * (t + 0.5) * dt));
      _frequency_transform_m(n, t) =
          dt * std::exp(-1i * 2.0 * constant::PI * (_frequencies(n) * t * dt));
    }
  }
}

void NffftFd::update() {
  // auto current_time_step{getCurrentTimeStep()};
  // updateXN(current_time_step);
  // updateXP(current_time_step);
  // updateYN(current_time_step);
  // updateYP(current_time_step);
  // updateZN(current_time_step);
  // updateZP(current_time_step);

  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nt{getCurrentTimeStep()};
  auto emf{getEMFInstance()};

  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto my_xn = -1 * emf->getEzFaceXCenter(li, j, k);
      auto mz_xn = emf->getEyFaceXCenter(li, j, k);
      auto jy_xn = emf->getHzFaceXCenter(li, j, k);
      auto jz_xn = -1 * emf->getHyFaceXCenter(li, j, k);

      auto my_xp = emf->getEzFaceXCenter(ri, j, k);
      auto mz_xp = -1 * emf->getEyFaceXCenter(ri, j, k);
      auto jy_xp = -1 * emf->getHzFaceXCenter(ri, j, k);
      auto jz_xp = emf->getHyFaceXCenter(ri, j, k);

      for (int n{0}; n < _number_frequencies; ++n) {
        _my_xn(n, j - lj, k - lk) += my_xn * _frequency_transform_m(n, nt);
        _mz_xn(n, j - lj, k - lk) += mz_xn * _frequency_transform_m(n, nt);
        _jy_xn(n, j - lj, k - lk) += jy_xn * _frequency_transform_j(n, nt);
        _jz_xn(n, j - lj, k - lk) += jz_xn * _frequency_transform_j(n, nt);

        _my_xp(n, j - lj, k - lk) += my_xp * _frequency_transform_m(n, nt);
        _mz_xp(n, j - lj, k - lk) += mz_xp * _frequency_transform_m(n, nt);
        _jy_xp(n, j - lj, k - lk) += jy_xp * _frequency_transform_j(n, nt);
        _jz_xp(n, j - lj, k - lk) += jz_xp * _frequency_transform_j(n, nt);
      }
    }
  }

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto mz_yn = -1 * emf->getExFaceYCenter(i, lj, k);
      auto mx_yn = emf->getEzFaceYCenter(i, lj, k);
      auto jz_yn = emf->getHxFaceYCenter(i, lj, k);
      auto jx_yn = -1 * emf->getHzFaceYCenter(i, lj, k);

      auto mz_yp = emf->getExFaceYCenter(i, rj, k);
      auto mx_yp = -1 * emf->getEzFaceYCenter(i, rj, k);
      auto jz_yp = -1 * emf->getHxFaceYCenter(i, rj, k);
      auto jx_yp = emf->getHzFaceYCenter(i, rj, k);

      for (int n{0}; n < _number_frequencies; ++n) {
        _mz_yn(n, k - lk, i - li) += mz_yn * _frequency_transform_m(n, nt);
        _mx_yn(n, k - lk, i - li) += mx_yn * _frequency_transform_m(n, nt);
        _jz_yn(n, k - lk, i - li) += jz_yn * _frequency_transform_j(n, nt);
        _jx_yn(n, k - lk, i - li) += jx_yn * _frequency_transform_j(n, nt);

        _mz_yp(n, k - lk, i - li) += mz_yp * _frequency_transform_m(n, nt);
        _mx_yp(n, k - lk, i - li) += mx_yp * _frequency_transform_m(n, nt);
        _jz_yp(n, k - lk, i - li) += jz_yp * _frequency_transform_j(n, nt);
        _jx_yp(n, k - lk, i - li) += jx_yp * _frequency_transform_j(n, nt);
      }
    }
  }

  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto mx_zn = -1 * emf->getEyFaceZCenter(i, j, lk);
      auto my_zn = emf->getExFaceZCenter(i, j, lk);
      auto jx_zn = emf->getHyFaceZCenter(i, j, lk);
      auto jy_zn = -1 * emf->getHxFaceZCenter(i, j, lk);

      auto mx_zp = emf->getEyFaceZCenter(i, j, rk);
      auto my_zp = -1 * emf->getExFaceZCenter(i, j, rk);
      auto jx_zp = -1 * emf->getHyFaceZCenter(i, j, rk);
      auto jy_zp = emf->getHxFaceZCenter(i, j, rk);

      for (int n{0}; n < _number_frequencies; ++n) {
        _mx_zn(n, i - li, j - lj) += mx_zn * _frequency_transform_m(n, nt);
        _my_zn(n, i - li, j - lj) += my_zn * _frequency_transform_m(n, nt);
        _jx_zn(n, i - li, j - lj) += jx_zn * _frequency_transform_j(n, nt);
        _jy_zn(n, i - li, j - lj) += jy_zn * _frequency_transform_j(n, nt);

        _mx_zp(n, i - li, j - lj) += mx_zp * _frequency_transform_m(n, nt);
        _my_zp(n, i - li, j - lj) += my_zp * _frequency_transform_m(n, nt);
        _jx_zp(n, i - li, j - lj) += jx_zp * _frequency_transform_j(n, nt);
        _jy_zp(n, i - li, j - lj) += jy_zp * _frequency_transform_j(n, nt);
      }
    }
  }
}

void NffftFd::updateXN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};

  for (auto j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto my{-0.5 * (getEz(li, j + 1, k) + getEz(li, j, k))};
      auto mz{0.5 * (getEy(li, j, k + 1) + getEy(li, j, k))};
      auto jy{0.25 * (getHz(li, j, k + 1) + getHz(li, j, k) +
                      getHz(li - 1, j, k + 1) + getHz(li - 1, j, k))};
      auto jz{-0.25 * (getHy(li, j + 1, k) + getHy(li, j, k) +
                       getHy(li - 1, j + 1, k) + getHy(li - 1, j, k))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _my_xn(n, j - lj, k - lk) +=
            my * _frequency_transform_m(n, current_time_step);
        _mz_xn(n, j - lj, k - lk) +=
            mz * _frequency_transform_m(n, current_time_step);
        _jy_xn(n, j - lj, k - lk) +=
            jy * _frequency_transform_j(n, current_time_step);
        _jz_xn(n, j - lj, k - lk) +=
            jz * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::updateXP(size_t current_time_step) {
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};

  for (auto j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto my{0.5 * (getEz(ri, j + 1, k) + getEz(ri, j, k))};
      auto mz{-0.5 * (getEy(ri, j, k + 1) + getEy(ri, j, k))};
      auto jy{-0.25 * (getHz(ri, j, k + 1) + getHz(ri, j, k) +
                       getHz(ri - 1, j, k + 1) + getHz(ri - 1, j, k))};
      auto jz{0.25 * (getHy(ri, j + 1, k) + getHy(ri, j, k) +
                      getHy(ri - 1, j + 1, k) + getHy(ri - 1, j, k))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _my_xp(n, j - lj, k - lk) +=
            my * _frequency_transform_m(n, current_time_step);
        _mz_xp(n, j - lj, k - lk) +=
            mz * _frequency_transform_m(n, current_time_step);
        _jy_xp(n, j - lj, k - lk) +=
            jy * _frequency_transform_j(n, current_time_step);
        _jz_xp(n, j - lj, k - lk) +=
            jz * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::updateYN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rk{getEndIndexZ()};

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto mz{-0.5 * (getEx(i, lj, k + 1) + getEx(i, lj, k))};
      auto mx{0.5 * (getEz(i + 1, lj, k) + getEz(i, lj, k))};
      auto jz{0.25 * (getHx(i + 1, lj, k) + getHx(i, lj, k) +
                      getHx(i + 1, lj - 1, k) + getHx(i, lj - 1, k))};
      auto jx{-0.25 * (getHz(i, lj, k + 1) + getHz(i, lj, k) +
                       getHz(i, lj - 1, k + 1) + getHz(i, lj - 1, k))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _mz_yn(n, k - lk, i - li) +=
            mz * _frequency_transform_m(n, current_time_step);
        _mx_yn(n, k - lk, i - li) +=
            mx * _frequency_transform_m(n, current_time_step);
        _jz_yn(n, k - lk, i - li) +=
            jz * _frequency_transform_j(n, current_time_step);
        _jx_yn(n, k - lk, i - li) +=
            jx * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::updateYP(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto mz{0.5 * (getEx(i, rj, k + 1) + getEx(i, rj, k))};
      auto mx{-0.5 * (getEz(i + 1, rj, k) + getEz(i, rj, k))};
      auto jz{-0.25 * (getHx(i + 1, rj, k) + getHx(i, rj, k) +
                       getHx(i + 1, rj - 1, k) + getHx(i, rj - 1, k))};
      auto jx{0.25 * (getHz(i, rj, k + 1) + getHz(i, rj, k) +
                      getHz(i, rj - 1, k + 1) + getHz(i, rj - 1, k))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _mz_yp(n, k - lk, i - li) +=
            mz * _frequency_transform_m(n, current_time_step);
        _mx_yp(n, k - lk, i - li) +=
            mx * _frequency_transform_m(n, current_time_step);
        _jz_yp(n, k - lk, i - li) +=
            jz * _frequency_transform_j(n, current_time_step);
        _jx_yp(n, k - lk, i - li) +=
            jx * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::updateZN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};

  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto mx{-0.5 * (getEy(i + 1, j, lk) + getEy(i, j, lk))};
      auto my{0.5 * (getEx(i, j + 1, lk) + getEx(i, j, lk))};
      auto jx{0.25 * (getHy(i, j + 1, lk) + getHy(i, j, lk) +
                      getHy(i, j + 1, lk - 1) + getHy(i, j, lk))};
      auto jy{-0.25 * (getHx(i + 1, j, lk) + getHx(i, j, lk) +
                       getHx(i + 1, j, lk - 1) + getHx(i, j, lk - 1))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _mx_zn(n, i - li, j - lj) +=
            mx * _frequency_transform_m(n, current_time_step);
        _my_zn(n, i - li, j - lj) +=
            my * _frequency_transform_m(n, current_time_step);
        _jx_zn(n, i - li, j - lj) +=
            jx * _frequency_transform_j(n, current_time_step);
        _jy_zn(n, i - li, j - lj) +=
            jy * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::updateZP(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};

  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto mx{0.5 * (getEy(i + 1, j, rk) + getEy(i, j, rk))};
      auto my{-0.5 * (getEx(i, j + 1, rk) + getEx(i, j, rk))};
      auto jx{-0.25 * (getHy(i, j + 1, rk) + getHy(i, j, rk) +
                       getHy(i, j + 1, rk - 1) + getHy(i, j, rk))};
      auto jy{0.25 * (getHx(i + 1, j, rk) + getHx(i, j, rk) +
                      getHx(i + 1, j, rk - 1) + getHx(i, j, rk - 1))};
      for (int n{0}; n < _number_frequencies; ++n) {
        _mx_zp(n, i - li, j - lj) +=
            mx * _frequency_transform_m(n, current_time_step);
        _my_zp(n, i - li, j - lj) +=
            my * _frequency_transform_m(n, current_time_step);
        _jx_zp(n, i - li, j - lj) +=
            jx * _frequency_transform_j(n, current_time_step);
        _jy_zp(n, i - li, j - lj) +=
            jy * _frequency_transform_j(n, current_time_step);
      }
    }
  }
}

void NffftFd::calculateFarfield() {
  for (int n{0}; n < _number_frequencies; ++n) {
    for (int t{0}; t < _number_theta; ++t) {
      for (int p{0}; p < _number_phi; ++p) {
        auto cos_t{std::cos(_theta(t))};
        auto sin_t{std::sin(_theta(t))};
        auto cos_p{std::cos(_phi(p))};
        auto sin_p{std::sin(_phi(p))};
        auto sin_t_cos_p{sin_t * cos_p};
        auto sin_t_sin_p{sin_t * sin_p};

        calculateFarfieldX(n, t, p, sin_t_cos_p, sin_t_sin_p, cos_t);
        calculateFarfieldY(n, t, p, sin_t_cos_p, sin_t_sin_p, cos_t);
        calculateFarfieldZ(n, t, p, sin_t_cos_p, sin_t_sin_p, cos_t);

        // xt::xtensor<double, 2> u{{cos_t * cos_p, cos_t * sin_p, -sin_t},
        //                          {-sin_p, cos_p, 0}};
        auto cos_t_cos_p{cos_t * cos_p};
        auto cos_t_sin_p{cos_t * sin_p};

        _a_theta(n, t, p) = cos_t_cos_p * _a_x(n, t, p) +
                            cos_t_sin_p * _a_y(n, t, p) - sin_t * _a_z(n, t, p);
        _a_phi(n, t, p) = -sin_p * _a_x(n, t, p) + cos_p * _a_y(n, t, p);
        _f_theta(n, t, p) = cos_t_cos_p * _f_x(n, t, p) +
                            cos_t_sin_p * _f_y(n, t, p) - sin_t * _f_z(n, t, p);
        _f_phi(n, t, p) = -sin_p * _f_x(n, t, p) + cos_p * _f_y(n, t, p);
      }
    }
  }
  using namespace std::complex_literals;
  auto coff{pow(_wave_number, 2) /
            (32 * constant::PI * constant::PI * constant::ETA_0)};
  _e_theta = (-1i * _wave_number) * (constant::ETA_0 * _a_theta + _f_phi) /
             (4 * constant::PI);
  _e_phi = (1i * _wave_number) * (-constant::ETA_0 * _a_phi + _f_theta) /
           (4 * constant::PI);
  _h_theta = (1i * _wave_number) * (_a_phi - _f_theta / constant::ETA_0) /
             (4 * constant::PI);
  _h_phi = (-1i * _wave_number) * (_a_theta + _f_phi / constant::ETA_0) /
           (4 * constant::PI);
  _power_theta = coff * pow(abs(constant::ETA_0 * _a_theta + _f_phi), 2);
  _power_phi = coff * pow(abs(-constant::ETA_0 * _a_phi + _f_theta), 2);
}

void NffftFd::calculateFarfieldX(int n, int t, int p, double sin_t_cos_p,
                                 double sin_t_sin_p, double cos_t) {
  using namespace std::complex_literals;
  const auto left_i{getStartIndexX()};  // don't use "li". to avoid mistakes
                                        // with imaginary number 1i
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const double dx_k{getDx() * _wave_number(n)};
  const double dy{getDy()};
  const double dz{getDz()};
  const auto dy_dz{dy * dz};  // area of rectangle

  std::complex<double> sum_jy{0, 0};
  std::complex<double> sum_jz{0, 0};
  std::complex<double> sum_my{0, 0};
  std::complex<double> sum_mz{0, 0};
  for (auto j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto delay_xn{dx_k *
                    ((left_i - 0) * sin_t_cos_p + (j - 0 + 0.5) * sin_t_sin_p +
                     (k - 0 + 0.5) * cos_t)};  // TODO(franzero)
      auto phase_delay_xn{std::exp(1i * delay_xn)};
      auto delay_xp{dx_k *
                    ((ri - 0) * sin_t_cos_p + (j - 0 + 0.5) * sin_t_sin_p +
                     (k - 0 + 0.5) * cos_t)};
      auto phase_delay_xp{std::exp(1i * delay_xp)};

      sum_jy += _jy_xn(n, j - lj, k - lk) * phase_delay_xn +
                _jy_xp(n, j - lj, k - lk) * phase_delay_xp;
      sum_jz += _jz_xn(n, j - lj, k - lk) * phase_delay_xn +
                _jz_xp(n, j - lj, k - lk) * phase_delay_xp;
      sum_my += _my_xn(n, j - lj, k - lk) * phase_delay_xn +
                _my_xp(n, j - lj, k - lk) * phase_delay_xp;
      sum_mz += _mz_xn(n, j - lj, k - lk) * phase_delay_xn +
                _mz_xp(n, j - lj, k - lk) * phase_delay_xp;
    }
  }
  _a_y(n, t, p) += sum_jy * dy_dz;
  _a_z(n, t, p) += sum_jz * dy_dz;
  _f_y(n, t, p) += sum_my * dy_dz;
  _f_z(n, t, p) += sum_mz * dy_dz;
}

void NffftFd::calculateFarfieldY(int n, int t, int p, double sin_t_cos_p,
                                 double sin_t_sin_p, double cos_t) {
  using namespace std::complex_literals;
  const auto left_i{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const double dx{getDx()};
  const double dy_k{getDy() * _wave_number(n)};
  const double dz{getDz()};
  const auto dz_dx{dz * dx};

  std::complex<double> sum_jz{0, 0};
  std::complex<double> sum_jx{0, 0};
  std::complex<double> sum_mz{0, 0};
  std::complex<double> sum_mx{0, 0};

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{left_i}; i < ri; ++i) {
      auto delay_yn{dy_k * ((i - 0 + 0.5) * sin_t_cos_p +
                            (lj - 0) * sin_t_sin_p + (k - 0 + 0.5) * cos_t)};
      auto phase_delay_yn{std::exp(1i * delay_yn)};
      auto delay_yp{dy_k * ((i - 0 + 0.5) * sin_t_cos_p +
                            (rj - 0) * sin_t_sin_p + (k - 0 + 0.5) * cos_t)};
      auto phase_delay_yp{std::exp(1i * delay_yp)};

      sum_jz += _jz_yn(n, k - lk, i - left_i) * phase_delay_yn +
                _jz_yp(n, k - lk, i - left_i) * phase_delay_yp;
      sum_jx += _jx_yn(n, k - lk, i - left_i) * phase_delay_yn +
                _jx_yp(n, k - lk, i - left_i) * phase_delay_yp;
      sum_mz += _mz_yn(n, k - lk, i - left_i) * phase_delay_yn +
                _mz_yp(n, k - lk, i - left_i) * phase_delay_yp;
      sum_mx += _mx_yn(n, k - lk, i - left_i) * phase_delay_yn +
                _mx_yp(n, k - lk, i - left_i) * phase_delay_yp;
    }
  }
  _a_z(n, t, p) += sum_jz * dz_dx;
  _a_x(n, t, p) += sum_jx * dz_dx;
  _f_z(n, t, p) += sum_mz * dz_dx;
  _f_x(n, t, p) += sum_mx * dz_dx;
}

void NffftFd::calculateFarfieldZ(int n, int t, int p, double sin_t_cos_p,
                                 double sin_t_sin_p, double cos_t) {
  using namespace std::complex_literals;
  const auto left_i{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const double dx{getDx()};
  const double dy{getDy()};
  const double dz_k{getDz() * _wave_number(n)};
  auto dx_dy{dx * dy};

  std::complex<double> sum_jx{0, 0};
  std::complex<double> sum_jy{0, 0};
  std::complex<double> sum_mx{0, 0};
  std::complex<double> sum_my{0, 0};

  for (auto j{lj}; j < rj; ++j) {
    for (auto i{left_i}; i < ri; ++i) {
      auto delay_zn{dz_k * ((i - 0 + 0.5) * sin_t_cos_p +
                            (j - 0 + 0.5) * sin_t_sin_p + (lk - 0) * cos_t)};
      auto phase_delay_zn{std::exp(1i * delay_zn)};
      auto delay_zp{dz_k * ((i - 0 + 0.5) * sin_t_cos_p +
                            (j - 0 + 0.5) * sin_t_sin_p + (rk - 0) * cos_t)};
      auto phase_delay_zp{std::exp(1i * delay_zp)};

      sum_jx += _jx_zn(n, i - left_i, j - lj) * phase_delay_zn +
                _jx_zp(n, i - left_i, j - lj) * phase_delay_zp;
      sum_jy += _jy_zn(n, i - left_i, j - lj) * phase_delay_zn +
                _jy_zp(n, i - left_i, j - lj) * phase_delay_zp;
      sum_mx += _mx_zn(n, i - left_i, j - lj) * phase_delay_zn +
                _mx_zp(n, i - left_i, j - lj) * phase_delay_zp;
      sum_my += _my_zn(n, i - left_i, j - lj) * phase_delay_zn +
                _my_zp(n, i - left_i, j - lj) * phase_delay_zp;
    }
  }

  _a_x(n, t, p) += sum_jx * dx_dy;
  _a_y(n, t, p) += sum_jy * dx_dy;
  _f_x(n, t, p) += sum_mx * dx_dy;
  _f_y(n, t, p) += sum_my * dx_dy;
}

void NffftFd::calculateFarfield1(
    const PointVector& r, SpatialIndex range_a_start, SpatialIndex range_a_end,
    SpatialIndex range_b_start, SpatialIndex range_b_end,
    SpatialIndex range_c_start, SpatialIndex range_c_end,
    const PointVector& offset) {
  using namespace std::complex_literals;
  auto sin_t_cos_p{r(0)};
  auto sin_t_sin_p{r(1)};
  auto cos_p{r(2)};
  for (auto i{range_a_start}; i < range_a_end; ++i) {
    for (auto j{range_b_start}; j < range_b_end; ++j) {
      for (auto k{range_c_start}; k < range_c_end; ++k) {
        auto delay{(i - 0 - offset(0) * sin_t_cos_p +
                    (j - 0 - offset(1)) * sin_t_sin_p +
                    (k - 0 - offset(2)) * cos_p)};
        auto phase_delay{std::exp(1i * delay)};
      }
    }
  }
}

void NffftFd::outputData() {
  calculateFarfield();
  const auto output_dir_path{getOutputDirPath()};
  if (!std::filesystem::exists(output_dir_path) ||
      !std::filesystem::is_directory(output_dir_path)) {
    try {
      std::filesystem::create_directory(output_dir_path);
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << output_dir_path
                << "\t Error:" << e.what() << '\n';
      return;
    }
  }
  // outputFarFieldParameters();

  xt::dump_npy(output_dir_path / "frequencies.npy", _frequencies);
  xt::dump_npy(output_dir_path / "theta.npy", _theta);
  xt::dump_npy(output_dir_path / "phi.npy", _phi);

  xt::dump_npy(output_dir_path / "e_theta.npy", _e_theta);
  xt::dump_npy(output_dir_path / "e_phi.npy", _e_phi);
  xt::dump_npy(output_dir_path / "h_theta.npy", _h_theta);
  xt::dump_npy(output_dir_path / "h_phi.npy", _h_phi);
  xt::dump_npy(output_dir_path / "n_theta.npy", _a_theta);
  xt::dump_npy(output_dir_path / "n_phi.npy", _a_phi);
  xt::dump_npy(output_dir_path / "l_theta.npy", _f_theta);
  xt::dump_npy(output_dir_path / "l_phi.npy", _f_phi);
  xt::dump_npy(output_dir_path / "power_theta.npy", _power_theta);
  xt::dump_npy(output_dir_path / "power_phi.npy", _power_phi);

  // for (int n{0}; n < _number_frequencies; ++n) {
  //   auto ghz{(_frequencies(n) / 1e9)};
  //   std::filesystem::path path{output_dir_path /
  //                              ("far_field_" + std::to_string(ghz) + "GHz")};
  //   try {
  //     std::filesystem::create_directory(path);
  //   } catch (std::exception e) {
  //     std::cerr << "Error: cannot create directory " << path
  //               << "\t Error:" << e.what() << '\n';
  //     return;
  //   }
  //   std::ofstream e_theta_data{path / ("e_theta.dat")};
  //   std::ofstream e_phi_data{path / ("e_phi_.dat")};
  //   std::ofstream h_theta_data{path / ("h_theta.dat")};
  //   std::ofstream h_phi_data{path / ("h_phi_.dat")};
  //   std::ofstream n_theta_data{path / ("n_theta.dat")};
  //   std::ofstream n_phi_data{path / ("n_phi_.dat")};
  //   std::ofstream l_theta_data{path / ("l_theta.dat")};
  //   std::ofstream l_phi_data{path / ("l_phi.dat")};
  //   std::ofstream power_theta_data{path / ("power_theta.dat")};
  //   std::ofstream power_phi_data{path / ("power_phi.dat")};

  //   for (int t{0}; t < _number_theta; ++t) {
  //     for (int p{0}; p < _number_phi; ++p) {
  //       e_theta_data << _e_theta(n, t, p) << "\t";
  //       e_phi_data << _e_phi(n, t, p) << "\t";
  //       h_theta_data << _h_theta(n, t, p) << "\t";
  //       h_phi_data << _h_phi(n, t, p) << "\t";
  //       n_theta_data << _a_theta(n, t, p) << "\t";
  //       n_phi_data << _a_phi(n, t, p) << "\t";
  //       l_theta_data << _f_theta(n, t, p) << "\t";
  //       l_phi_data << _f_phi(n, t, p) << "\t";

  //       power_theta_data << _power_theta(n, t, p) << "\t";
  //       power_phi_data << _power_phi(n, t, p) << "\t";
  //     }
  //     e_theta_data << std::endl;
  //     e_phi_data << std::endl;
  //     h_theta_data << std::endl;
  //     h_phi_data << std::endl;
  //     n_theta_data << std::endl;
  //     n_phi_data << std::endl;
  //     l_theta_data << std::endl;
  //     l_phi_data << std::endl;
  //     power_theta_data << std::endl;
  //     power_phi_data << std::endl;
  //   }
  //   e_theta_data.close();
  //   e_phi_data.close();
  //   h_theta_data.close();
  //   h_phi_data.close();
  //   power_theta_data.close();
  //   power_phi_data.close();
  // }
}

void NffftFd::outputFarFieldParameters() {
  // std::ofstream far_field_parameter_writer{getOutputDirPath() /
  //                                          "far_field_parameter.dat"};
  // far_field_parameter_writer << "number of frequencies: " <<
  // _number_frequencies
  //                            << "\n";
  // for (const auto& e : _frequencies) {
  //   far_field_parameter_writer << e << "\t";
  // }
  // far_field_parameter_writer << std::endl;
  // far_field_parameter_writer << "number of theta: " << _number_theta << "\n";
  // for (const auto& e : _theta) {
  //   far_field_parameter_writer << e << "\t";
  // }
  // far_field_parameter_writer << std::endl;
  // far_field_parameter_writer << "number of phi: " << _number_phi << "\n";
  // for (const auto& e : _phi) {
  //   far_field_parameter_writer << e << "\t";
  // }
  // far_field_parameter_writer << "\nthe intrinsic impedance: "
  //                            << constant::ETA_0;
  // far_field_parameter_writer.close();
}
}  // namespace xfdtd