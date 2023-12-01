#include "nffft/nffft_broadband.h"

#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xnpy.hpp>

#include "grid/grid_space.h"
#include "nffft/nffft.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
NffftBroadBand::NffftBroadBand(SpatialIndex distance_x, SpatialIndex distance_y,
                               SpatialIndex distance_z, double theta,
                               double phi, std::string output_dir_path)
    : NFFFT(distance_x, distance_y, distance_z, std::move(output_dir_path)),
      _theta{theta},
      _phi{phi},
      _fairfield_vector{sin(theta) * cos(phi), sin(theta) * sin(phi),
                        cos(theta)} {}

void NffftBroadBand::init(std::shared_ptr<const GridSpace> grid_space,
                          std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                          std::shared_ptr<const EMF> emf) {
  defaultInit(std::move(grid_space), std::move(fdtd_basic_coff),
              std::move(emf));

  _number_samples = getTotalTimeSteps() * 2 - 1;  // make sure it is odd number
  _sample_rate = 1 / getDt();

  _jy_xn.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jy_xp.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xn.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _jz_xp.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _jz_yn.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _jz_yp.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _jx_yn.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _jx_yp.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});

  _jx_zn.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _jx_zp.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _jy_zn.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _jy_zp.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});

  _my_xn.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _my_xp.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xn.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});
  _mz_xp.resize({_number_samples, 1, static_cast<size_t>(getNy()),
                 static_cast<size_t>(getNz())});

  _mz_yn.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _mz_yp.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _mx_yn.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});
  _mx_yp.resize({_number_samples, static_cast<size_t>(getNx()), 1,
                 static_cast<size_t>(getNz())});

  _mx_zn.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _mx_zp.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _my_zn.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
  _my_zp.resize({_number_samples, static_cast<size_t>(getNx()),
                 static_cast<size_t>(getNy()), 1});
}

void NffftBroadBand::update() {
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
      _my_xn(nt, 0, j - lj, k - lk) = -1.0 * emf->getEzFaceXCenter(li, j, k);
      _mz_xn(nt, 0, j - lj, k - lk) = emf->getEyFaceXCenter(li, j, k);
      _jy_xn(nt, 0, j - lj, k - lk) = emf->getHzFaceXCenter(li, j, k);
      _jz_xn(nt, 0, j - lj, k - lk) = -1.0 * emf->getHyFaceXCenter(li, j, k);

      _my_xp(nt, 0, j - lj, k - lk) = emf->getEzFaceXCenter(ri, j, k);
      _mz_xp(nt, 0, j - lj, k - lk) = -1.0 * emf->getEyFaceXCenter(ri, j, k);
      _jy_xp(nt, 0, j - lj, k - lk) = -1.0 * emf->getHzFaceXCenter(ri, j, k);
      _jz_xp(nt, 0, j - lj, k - lk) = emf->getHyFaceXCenter(ri, j, k);
    }
  }

  // Y
  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      _mz_yn(nt, i - li, 0, k - lk) = -1.0 * emf->getExFaceYCenter(i, lj, k);
      _mx_yn(nt, i - li, 0, k - lk) = emf->getEzFaceYCenter(i, lj, k);
      _jz_yn(nt, i - li, 0, k - lk) = emf->getHxFaceYCenter(i, lj, k);
      _jx_yn(nt, i - li, 0, k - lk) = -1.0 * emf->getHzFaceYCenter(i, lj, k);

      _mz_yp(nt, i - li, 0, k - lk) = emf->getExFaceYCenter(i, rj, k);
      _mx_yp(nt, i - li, 0, k - lk) = -1.0 * emf->getEzFaceYCenter(i, rj, k);
      _jz_yp(nt, i - li, 0, k - lk) = -1.0 * emf->getHxFaceYCenter(i, rj, k);
      _jx_yp(nt, i - li, 0, k - lk) = emf->getHzFaceYCenter(i, rj, k);
    }
  }

  // Z
  for (SpatialIndex i{li}; i < ri; ++i) {
    for (SpatialIndex j{lj}; j < rj; ++j) {
      _mx_zn(nt, i - li, j - lj, 0) = -1.0 * emf->getEyFaceZCenter(i, j, lk);
      _my_zn(nt, i - li, j - lj, 0) = emf->getExFaceZCenter(i, j, lk);
      _jx_zn(nt, i - li, j - lj, 0) = emf->getHyFaceZCenter(i, j, lk);
      _jy_zn(nt, i - li, j - lj, 0) = -1.0 * emf->getHxFaceZCenter(i, j, lk);

      _mx_zp(nt, i - li, j - lj, 0) = emf->getEyFaceZCenter(i, j, rk);
      _my_zp(nt, i - li, j - lj, 0) = -1.0 * emf->getExFaceZCenter(i, j, rk);
      _jx_zp(nt, i - li, j - lj, 0) = -1.0 * emf->getHyFaceZCenter(i, j, rk);
      _jy_zp(nt, i - li, j - lj, 0) = emf->getHxFaceZCenter(i, j, rk);
    }
  }
}

void NffftBroadBand::outputData() {
  calculateFarField();
  auto output_path{std::filesystem::path{getOutputDirPath()}};
  if (!std::filesystem::exists(output_path) &&
      !std::filesystem::is_directory(output_path)) {
    try {
      std::filesystem::create_directories(output_path);
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << output_path.string()
                << '\n';
      return;
    }
  }

  outputFarFieldParameters();
  using namespace std::complex_literals;
  if (getGridSpace()->getDimension() == GridSpace::Dimension::TWO_DIMENSION) {
    // Debug
    auto far_ez{0.5 * xt::sqrt(1i * _wave_number) *
                (_l_phi + constant::ETA_0 * _n_theta)};
    xt::dump_npy(std::filesystem::path{output_path / "EZ.npy"}.string(),
                 far_ez);
    return;
  }

  // xt::dump_npy(output_dir_path / "e_theta.npy", _e_theta);
  // xt::dump_npy(output_dir_path / "e_phi.npy", _e_phi);
  // xt::dump_npy(output_dir_path / "h_theta.npy", _h_theta);
  // xt::dump_npy(output_dir_path / "h_phi.npy", _h_phi);
  // xt::dump_npy(output_dir_path / "power_theta.npy", _power_theta);
  // xt::dump_npy(output_dir_path / "power_phi.npy", _power_phi);

  xt::dump_npy(std::filesystem::path{output_path / "e_theta.npy"}.string(),
               _e_theta);
  xt::dump_npy(std::filesystem::path{output_path / "e_phi.npy"}.string(),
               _e_phi);
  xt::dump_npy(std::filesystem::path{output_path / "h_theta.npy"}.string(),
               _h_theta);
  xt::dump_npy(std::filesystem::path{output_path / "h_phi.npy"}.string(),
               _h_phi);
  xt::dump_npy(std::filesystem::path{output_path / "power_theta.npy"}.string(),
               _power_theta);
  xt::dump_npy(std::filesystem::path{output_path / "power_phi.npy"}.string(),
               _power_phi);
}

void NffftBroadBand::updateXN(size_t current_time_step) {
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
      _my_xn(nt, 0, j - lj, k - lk) = my;
      _mz_xn(nt, 0, j - lj, k - lk) = mz;
      _jy_xn(nt, 0, j - lj, k - lk) = jy;
      _jz_xn(nt, 0, j - lj, k - lk) = jz;
    }
  }
}

void NffftBroadBand::updateXP(size_t current_time_step) {
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
      _my_xp(nt, 0, j - lj, k - lk) = my;
      _mz_xp(nt, 0, j - lj, k - lk) = mz;
      _jy_xp(nt, 0, j - lj, k - lk) = jy;
      _jz_xp(nt, 0, j - lj, k - lk) = jz;
    }
  }
}

void NffftBroadBand::updateYN(size_t current_time_step) {
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rk{getEndIndexZ()};
  const auto nt{current_time_step};

  for (SpatialIndex k{lk}; k < rk; ++k) {
    for (SpatialIndex i{li}; i < ri; ++i) {
      auto mz{-0.5 * (getEx(i, lj, k + 1) + getEx(i, lj, k))};
      auto mx{0.5 * (getEz(i + 1, lj, k) + getEz(i, lj, k))};
      auto jz{0.25 * (getHx(i + 1, lj, k) + getHx(i, lj, k) +
                      getHx(i + 1, lj - 1, k) + getHx(i, lj - 1, k))};
      auto jx{-0.25 * (getHz(i, lj, k + 1) + getHz(i, lj, k) +
                       getHz(i, lj - 1, k + 1) + getHz(i, lj - 1, k))};
      _mz_yn(nt, i - li, 0, k - lk) = mz;
      _mx_yn(nt, i - li, 0, k - lk) = mx;
      _jz_yn(nt, i - li, 0, k - lk) = jz;
      _jx_yn(nt, i - li, 0, k - lk) = jx;
    }
  }
}

void NffftBroadBand::updateYP(size_t current_time_step) {
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
      _mz_yp(nt, i - li, 0, k - lk) = mz;
      _mx_yp(nt, i - li, 0, k - lk) = mx;
      _jz_yp(nt, i - li, 0, k - lk) = jz;
      _jx_yp(nt, i - li, 0, k - lk) = jx;
    }
  }
}

void NffftBroadBand::updateZN(size_t current_time_step) {
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
                      getHy(i, j + 1, lk - 1) + getHy(i, j, lk - 1))};
      auto jy{-0.25 * (getHx(i + 1, j, lk) + getHx(i, j, lk) +
                       getHx(i + 1, j, lk - 1) + getHx(i, j, lk - 1))};
      _mx_zn(nt, i - li, j - lj, 0) = mx;
      _my_zn(nt, i - li, j - lj, 0) = my;
      _jx_zn(nt, i - li, j - lj, 0) = jx;
      _jy_zn(nt, i - li, j - lj, 0) = jy;
    }
  }
}

void NffftBroadBand::updateZP(size_t current_time_step) {
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
                       getHy(i, j + 1, rk - 1) + getHy(i, j, rk - 1))};
      auto jy{0.25 * (getHx(i + 1, j, rk) + getHx(i, j, rk) +
                      getHx(i + 1, j, rk - 1) + getHx(i, j, rk - 1))};

      _mx_zp(nt, i - li, j - lj, 0) = mx;
      _my_zp(nt, i - li, j - lj, 0) = my;
      _jx_zp(nt, i - li, j - lj, 0) = jx;
      _jy_zp(nt, i - li, j - lj, 0) = jy;
    }
  }
}

void NffftBroadBand::calculateFarField() {
  std::cout << "Start calculating far field" << '\n';
  const auto fs{1 / getDt()};
  _frequencies = xt::linspace(-fs / 2, fs / 2, _number_samples);
  _wave_number = 2 * constant::PI * _frequencies / constant::C_0;

  calculateLNXYZ();

  auto cos_t{std::cos(_theta)};
  auto sin_t{std::sin(_theta)};
  auto cos_p{std::cos(_phi)};
  auto sin_p{std::sin(_phi)};
  auto sin_t_cos_p{sin_t * cos_p};
  auto sin_t_sin_p{sin_t * sin_p};

  _n_theta = cos_t * cos_p * _n_x + cos_t * sin_p * _n_y - sin_t * _n_z;
  _n_phi = -sin_p * _n_x + cos_p * _n_y;
  _l_theta = cos_t * cos_p * _l_x + cos_t * sin_p * _l_y - sin_t * _l_z;
  _l_phi = -sin_p * _l_x + cos_p * _l_y;
  using namespace std::complex_literals;
  _e_theta = (-1i * _wave_number) * (constant::ETA_0 * _n_theta + _l_phi) /
             (4 * constant::PI);
  _e_phi = (1i * _wave_number) * (-constant::ETA_0 * _n_phi + _l_theta) /
           (4 * constant::PI);
  _h_theta = (1i * _wave_number) * (_n_phi - _l_theta / constant::ETA_0) /
             (4 * constant::PI);
  _h_phi = (-1i * _wave_number) * (_n_theta + _l_phi / constant::ETA_0) /
           (4 * constant::PI);
  auto coff{pow(_wave_number, 2) /
            (32 * constant::PI * constant::PI * constant::ETA_0)};
  _power_theta = coff * pow(abs(constant::ETA_0 * _n_theta + _l_phi), 2);
  _power_phi = coff * pow(abs(-constant::ETA_0 * _n_phi + _l_theta), 2);
  std::cout << "Finish calculating far field" << '\n';
}

void NffftBroadBand::calculateLNXYZ() {
  const auto dx{getDx()};
  const auto dy{getDy()};
  const auto dz{getDz()};
  const auto li{getStartIndexX()};
  const auto lj{getStartIndexY()};
  const auto lk{getStartIndexZ()};
  const auto ri{getEndIndexX()};
  const auto rj{getEndIndexY()};
  const auto rk{getEndIndexZ()};
  const auto nx{getNx()};
  const auto ny{getNy()};
  const auto nz{getNz()};

  auto calculate_delay = [&](SpatialIndex li, SpatialIndex ri, SpatialIndex lj,
                             SpatialIndex rj, SpatialIndex lk, SpatialIndex rk,
                             double offset_x, double offset_y,
                             double offset_z) -> xt::xarray<double> {
    if (li == ri) {
      ri += 1;
    }
    if (lj == rj) {
      rj += 1;
    }
    if (lk == rk) {
      rk += 1;
    }
    auto i_range{xt::arange<double>(li, ri, 1)};
    auto j_range{xt::arange<double>(lj, rj, 1)};
    auto k_range{xt::arange<double>(lk, rk, 1)};
    auto [x, y, z]{xt::meshgrid(i_range + offset_x - getCenterIndexX(),
                                j_range + offset_y - getCenterIndexY(),
                                k_range + offset_z - getCenterIndexZ())};
    auto x_r{x * _fairfield_vector(0) * dx};
    auto y_r{y * _fairfield_vector(1) * dy};
    auto z_r{z * _fairfield_vector(2) * dz};

    return {(x_r + y_r + z_r)};
  };

  auto delay_xn{calculate_delay(li, li, lj, rj, lk, rk, 0, 0.5, 0.5)};
  auto delay_xp{calculate_delay(ri, ri, lj, rj, lk, rk, 0, 0.5, 0.5)};
  auto delay_yn{calculate_delay(li, ri, lj, lj, lk, rk, 0.5, 0, 0.5)};
  auto delay_yp{calculate_delay(li, ri, rj, rj, lk, rk, 0.5, 0, 0.5)};
  auto delay_zn{calculate_delay(li, ri, lj, rj, lk, lk, 0.5, 0.5, 0)};
  auto delay_zp{calculate_delay(li, ri, lj, rj, rk, rk, 0.5, 0.5, 0)};

  calculateAuxiliary(_n_y, _n_z, _l_y, _l_z, _jy_xn, _jz_xn, _my_xn, _mz_xn, 1,
                     ny, nz, delay_xn, dy * dz);
  calculateAuxiliary(_n_y, _n_z, _l_y, _l_z, _jy_xp, _jz_xp, _my_xp, _mz_xp, 1,
                     ny, nz, delay_xp, dy * dz);
  calculateAuxiliary(_n_x, _n_z, _l_x, _l_z, _jx_yn, _jz_yn, _mx_yn, _mz_yn, nx,
                     1, nz, delay_yn, dx * dz);
  calculateAuxiliary(_n_x, _n_z, _l_x, _l_z, _jx_yp, _jz_yp, _mx_yp, _mz_yp, nx,
                     1, nz, delay_yp, dx * dz);
  if (nz <= 1) {
    return;
  }
  calculateAuxiliary(_n_x, _n_y, _l_x, _l_y, _jx_zn, _jy_zn, _mx_zn, _my_zn, nx,
                     ny, 1, delay_zn, dx * dy);
  calculateAuxiliary(_n_x, _n_y, _l_x, _l_y, _jx_zp, _jy_zp, _mx_zp, _my_zp, nx,
                     ny, 1, delay_zp, dx * dy);
}

void NffftBroadBand::calculateAuxiliary(EFFA &n_a, EFFA &n_b, EFFA &l_a,
                                        EFFA &l_b, const EFTA &j_a,
                                        const EFTA &j_b, EFTA &m_a, EFTA &m_b,
                                        int x, int y, int z,
                                        const xt::xarray<double> &delay,
                                        double ds) {
  using namespace std::complex_literals;
  for (int i{0}; i < x; ++i) {
    for (int j{0}; j < y; ++j) {
      for (int k{0}; k < z; ++k) {
        auto phase_shift{xt::exp(1i * delay(i, j, k) * _wave_number)};
        n_a += ds * phase_shift *
               xt::fftw::fftshift(
                   xt::fftw::fft(xt::view(j_a, xt::all(), i, j, k)));
        n_b += ds * phase_shift *
               xt::fftw::fftshift(
                   xt::fftw::fft(xt::view(j_b, xt::all(), i, j, k)));
        l_a += ds * phase_shift *
               xt::fftw::fftshift(
                   xt::fftw::fft(xt::view(m_a, xt::all(), i, j, k)));
        l_b += ds * phase_shift *
               xt::fftw::fftshift(
                   xt::fftw::fft(xt::view(m_b, xt::all(), i, j, k)));
      }
    }
  }
}

void NffftBroadBand::outputFarFieldParameters() {
  std::ofstream far_field_parameter_writer{std::filesystem::path{
      std::filesystem::path{getOutputDirPath()} / "far_field_parameter.dat"}
                                               .string()};
  far_field_parameter_writer << "Fs:\t" << 1 / getDt() << '\n';
  far_field_parameter_writer << "number of samples:\t" << _number_samples
                             << '\n';
  far_field_parameter_writer << "Theta:\n" << _theta << '\n';
  far_field_parameter_writer << "Phi:\n" << _phi << '\n';
  far_field_parameter_writer.close();
}

}  // namespace xfdtd
