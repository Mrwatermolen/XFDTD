#include "nffft/nffft_2d_test.h"

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include "nffft/nffft.h"

namespace xfdtd {
NFFFT2DTEST::NFFFT2DTEST(SpatialIndex distance_x, SpatialIndex distance_y,
                         SpatialIndex direction_z, double far_theta,
                         double far_phi, std::string output_dir_path)
    : NFFFT(distance_x, distance_y, direction_z, std::move(output_dir_path)),
      _farfield_theta{far_theta},
      _farfield_phi{far_phi} {}

void NFFFT2DTEST::init(std::shared_ptr<const GridSpace> grid_space,
                       std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                       std::shared_ptr<const EMF> emf) {
  defaultInit(std::move(grid_space), std::move(fdtd_basic_coff),
              std::move(emf));

  _jy_xn.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jy_xp.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xn.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xp.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _jz_yn.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jz_yp.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jx_yn.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _jx_yp.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});

  _jx_zn.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jx_zp.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jy_zn.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _jy_zp.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});

  _my_xn.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _my_xp.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xn.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xp.resize({getTotalTimeSteps(), static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _mz_yn.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mz_yp.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mx_yn.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});
  _mx_yp.resize({getTotalTimeSteps(), static_cast<size_t>(getNz()),
                 static_cast<size_t>(getNx())});

  _mx_zn.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _mx_zp.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _my_zn.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
  _my_zp.resize({getTotalTimeSteps(), static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy())});
}

void NFFFT2DTEST::update() {
  auto current_time_step{getCurrentTimeStep()};
  updateXN(current_time_step);
  updateXP(current_time_step);
  updateYN(current_time_step);
  updateYP(current_time_step);
}

void NFFFT2DTEST::outputData() { calculateFarfield(); }

void NFFFT2DTEST::updateXN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto my{-0.5 * (getEz(li, j + 1, k) + getEz(li, j, k))};
      auto mz{0.5 * (getEy(li, j, k + 1) + getEy(li, j, k))};
      auto jy{0.25 * (getHz(li, j, k + 1) + getHz(li, j, k) +
                      getHz(li - 1, j, k + 1) + getHz(li - 1, j, k))};
      auto jz{-0.25 * (getHy(li, j + 1, k) + getHy(li, j, k) +
                       getHy(li - 1, j + 1, k) + getHy(li - 1, j, k))};
      _my_xn(nt, j - lj, k - lk) = my;
      _mz_xn(nt, j - lj, k - lk) = mz;
      _jy_xn(nt, j - lj, k - lk) = jy;
      _jz_xn(nt, j - lj, k - lk) = jz;
    }
  }
}

void NFFFT2DTEST::updateXP(size_t current_time_step) {
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (SpatialIndex j{lj}; j < rj; ++j) {
    for (auto k{lk}; k < rk; ++k) {
      auto my{0.5 * (getEz(ri, j + 1, k) + getEz(ri, j, k))};
      auto mz{-0.5 * (getEy(ri, j, k + 1) + getEy(ri, j, k))};
      auto jy{-0.25 * (getHz(ri, j, k + 1) + getHz(ri, j, k) +
                       getHz(ri - 1, j, k + 1) + getHz(ri - 1, j, k))};
      auto jz{0.25 * (getHy(ri, j + 1, k) + getHy(ri, j, k) +
                      getHy(ri - 1, j + 1, k) + getHy(ri - 1, j, k))};
      _my_xp(nt, j - lj, k - lk) = my;
      _mz_xp(nt, j - lj, k - lk) = mz;
      _jy_xp(nt, j - lj, k - lk) = jy;
      _jz_xp(nt, j - lj, k - lk) = jz;
    }
  }
}

void NFFFT2DTEST::updateYN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto mz{-0.5 * (getEx(i, lj, k + 1) + getEx(i, lj, k))};
      auto mx{0.5 * (getEz(i + 1, lj, k) + getEz(i, lj, k))};
      auto jz{0.25 * (getHx(i + 1, lj, k) + getHx(i, lj, k) +
                      getHx(i + 1, lj - 1, k) + getHx(i, lj - 1, k))};
      auto jx{-0.25 * (getHz(i, lj, k + 1) + getHz(i, lj, k) +
                       getHz(i, lj - 1, k + 1) + getHz(i, lj - 1, k))};
      _mz_yn(nt, k - lk, i - li) = mz;
      _mx_yn(nt, k - lk, i - li) = mx;
      _jz_yn(nt, k - lk, i - li) = jz;
      _jx_yn(nt, k - lk, i - li) = jx;
    }
  }
}

void NFFFT2DTEST::updateYP(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (auto k{lk}; k < rk; ++k) {
    for (auto i{li}; i < ri; ++i) {
      auto mz{0.5 * (getEx(i, rj, k + 1) + getEx(i, rj, k))};
      auto mx{-0.5 * (getEz(i + 1, rj, k) + getEz(i, rj, k))};
      auto jz{-0.25 * (getHx(i + 1, rj, k) + getHx(i, rj, k) +
                       getHx(i + 1, rj - 1, k) + getHx(i, rj - 1, k))};
      auto jx{0.25 * (getHz(i, rj, k + 1) + getHz(i, rj, k) +
                      getHz(i, rj - 1, k + 1) + getHz(i, rj - 1, k))};
      _mz_yp(nt, k - lk, i - li) = mz;
      _mx_yp(nt, k - lk, i - li) = mx;
      _jz_yp(nt, k - lk, i - li) = jz;
      _jx_yp(nt, k - lk, i - li) = jx;
    }
  }
}

void NFFFT2DTEST::updateZN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto nt{current_time_step};

  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto mx{-0.5 * (getEy(i + 1, j, lk) + getEy(i, j, lk))};
      auto my{0.5 * (getEx(i, j + 1, lk) + getEx(i, j, lk))};
      auto jx{0.25 * (getHy(i, j + 1, lk) + getHy(i, j, lk) +
                      getHy(i, j + 1, lk - 1) + getHy(i, j, lk))};
      auto jy{-0.25 * (getHx(i + 1, j, lk) + getHx(i, j, lk) +
                       getHx(i + 1, j, lk - 1) + getHx(i, j, lk - 1))};
      _mx_zn(nt, i - li, j - lj) = mx;
      _my_zn(nt, i - li, j - lj) = my;
      _jx_zn(nt, i - li, j - lj) = jx;
      _jy_zn(nt, i - li, j - lj) = jy;
    }
  }
}

void NFFFT2DTEST::updateZP(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (auto i{li}; i < ri; ++i) {
    for (auto j{lj}; j < rj; ++j) {
      auto mx{0.5 * (getEy(i + 1, j, rk) + getEy(i, j, rk))};
      auto my{-0.5 * (getEx(i, j + 1, rk) + getEx(i, j, rk))};
      auto jx{-0.25 * (getHy(i, j + 1, rk) + getHy(i, j, rk) +
                       getHy(i, j + 1, rk - 1) + getHy(i, j, rk))};
      auto jy{0.25 * (getHx(i + 1, j, rk) + getHx(i, j, rk) +
                      getHx(i + 1, j, rk - 1) + getHx(i, j, rk - 1))};
      _mx_zp(nt, i - li, j - lj) = mx;
      _my_zp(nt, i - li, j - lj) = my;
      _jx_zp(nt, i - li, j - lj) = jx;
      _jy_zp(nt, i - li, j - lj) = jy;
    }
  }
}

void NFFFT2DTEST::calculateFarfield() {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto dx{getDx()};
  const auto dy{getDy()};
  const auto dz{getDz()};
  const auto dt{getDt()};
  // const auto center_x{getOutputBoundaryCenterX()};
  // const auto center_y{getOutputBoundaryCenterY()};
  // const auto center_z{getOutputBoundaryCenterZ()};
  const auto center_x{0};
  const auto center_y{0};
  const auto center_z{0};
  const auto phi{getFarfieldPhi()};
  const auto cos_phi{cos(phi)};
  const auto sin_phi{sin(phi)};
  const auto dy_dz_sinth{dy * dz};
  const auto dy_dz_sinphi{dy * dz * sin_phi};
  const auto dy_dz_cosphi{dy * dz * cos_phi};
  const auto dz_dx_sinth{dz * dx};
  const auto dz_dx_sinphi{dz * dx * sin_phi};
  const auto dz_dx_cosphi{dz * dx * cos_phi};
  auto num_frequency{getTotalTimeSteps() * 2 - 1};
  // TODO(franzero): refactor this
  if (_jz_xn.shape(0) < num_frequency) {
    auto size{num_frequency - _jz_xn.shape(0)};
    auto additive_shape{_jz_xn.shape()};
    additive_shape[0] = size;
    _jz_xn = xt::concatenate(
        xt::xtuple(_jz_xn, xt::zeros<double>(additive_shape)), 0);
    _jz_xp = xt::concatenate(
        xt::xtuple(_jz_xp, xt::zeros<double>(additive_shape)), 0);
    _my_xn = xt::concatenate(
        xt::xtuple(_my_xn, xt::zeros<double>(additive_shape)), 0);
    _my_xp = xt::concatenate(
        xt::xtuple(_my_xp, xt::zeros<double>(additive_shape)), 0);
    additive_shape = _jz_yn.shape();
    additive_shape[0] = size;
    _jz_yn = xt::concatenate(
        xt::xtuple(_jz_yn, xt::zeros<double>(additive_shape)), 0);
    _jz_yp = xt::concatenate(
        xt::xtuple(_jz_yp, xt::zeros<double>(additive_shape)), 0);
    _mx_yn = xt::concatenate(
        xt::xtuple(_mx_yn, xt::zeros<double>(additive_shape)), 0);
    _mx_yp = xt::concatenate(
        xt::xtuple(_mx_yp, xt::zeros<double>(additive_shape)), 0);
  }
  auto fs{1 / dt};
  double step_frequency{1e8};
  auto frequencies{xt::linspace(-fs / 2, fs / 2, num_frequency)};
  auto k_wave{2 * constant::PI * frequencies / constant::C_0};

  auto fixed_delay_xn{(dx * (li - center_x) * cos_phi)};
  auto fixed_delay_xp{(dx * (ri - center_x) * cos_phi)};

  auto fixed_delay_yn{(dy * (lj - center_y) * sin_phi)};
  auto fixed_delay_yp{(dy * (rj - center_y) * sin_phi)};

  xt::xarray<std::complex<double>> n_theta{std::vector<size_t>{num_frequency}};
  xt::xarray<std::complex<double>> l_phi{std::vector<size_t>{num_frequency}};
  std::complex<double> i1{0, 1};
  SpatialIndex k{0};

  for (SpatialIndex i{li}; i < ri; ++i) {
    auto delay_x{dx * (i - center_x + 0.5) * cos_phi};
    auto exp_jk_rpr_yn{xt::exp(i1 * k_wave * (delay_x + fixed_delay_yn))};
    auto exp_jk_rpr_yp{xt::exp(i1 * k_wave * (delay_x + fixed_delay_yp))};

    xt::xarray<double> temp{xt::view(_jz_yn, xt::all(), k, i - li)};
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinth * exp_jk_rpr_yn;

    temp = xt::view(_jz_yp, xt::all(), k, i - li);
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinth * exp_jk_rpr_yp;

    temp = xt::view(_mx_yn, xt::all(), k, i - li);
    l_phi +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinphi * exp_jk_rpr_yn;

    temp = xt::view(_mx_yp, xt::all(), k, i - li);
    l_phi +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinphi * exp_jk_rpr_yp;
  }

  for (SpatialIndex j{lj}; j < rj; ++j) {
    // direction: XN XP
    auto delay_y{dy * (j - center_y + 0.5) * sin_phi};
    auto exp_jk_rpr_xn{xt::exp(i1 * k_wave * (delay_y + fixed_delay_xn))};
    auto exp_jk_rpr_xp{xt::exp(i1 * k_wave * (delay_y + fixed_delay_xp))};

    xt::xarray<double> temp{xt::view(_jz_xn, xt::all(), j - lj, k)};

    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_sinth * exp_jk_rpr_xn;

    temp = xt::view(_jz_xp, xt::all(), j - lj, k);
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_sinth * exp_jk_rpr_xp;

    temp = xt::view(_my_xn, xt::all(), j - lj, k);
    l_phi +=
        xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_cosphi * exp_jk_rpr_xn;

    temp = xt::view(_my_xp, xt::all(), j - lj, k);
    l_phi +=
        xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_cosphi * exp_jk_rpr_xp;
  }

  auto far_ez{0.5 * xt::sqrt(i1 * k_wave) *
              (l_phi + constant::ETA_0 * n_theta)};
  auto mangitude{xt::abs(far_ez)};
  std::ofstream fout;
  fout.open("visualizing_data/far_ez.dat");
  for (auto&& e : mangitude) {
    fout << e << " ";
  }
  fout.close();
}
}  // namespace xfdtd