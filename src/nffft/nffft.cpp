#include "nffft/nffft.h"

#include <complex>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"

namespace xfdtd {

NFFFT::NFFFT(SpatialIndex distance_x, SpatialIndex distance_y,
             SpatialIndex direction_z, double far_tehta, double far_phi)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{direction_z},
      _farfield_theta{far_tehta},
      _farfield_phi{far_phi} {}

void NFFFT::init(std::unique_ptr<GridBox> output_box, std::shared_ptr<EMF> emf,
                 size_t total_time_steps, double dt, double dx, double dy,
                 double dz) {
  if (output_box == nullptr) {
    throw std::runtime_error("Output box instance is not set.");
  }
  if (emf == nullptr) {
    throw std::runtime_error("EMF instance is not set.");
  }

  _output_box = std::move(output_box);
  _emf = std::move(emf);
  _total_time_steps = total_time_steps;
  _dt = dt;
  _dx = dx;
  _dy = dy;
  _dz = dz;

  _jez_xn = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNy(), getOutputBoundaryNz()));
  _jmy_xn = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNy(), getOutputBoundaryNz()));

  _jez_xp = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNy(), getOutputBoundaryNz()));
  _jmy_xp = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNy(), getOutputBoundaryNz()));

  _jez_yn = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNz(), getOutputBoundaryNx()));
  _jmx_yn = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNz(), getOutputBoundaryNx()));

  _jez_yp = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNz(), getOutputBoundaryNx()));
  _jmx_yp = std::move(allocateDoubleArray3D(
      getTotalTimeSteps(), getOutputBoundaryNz(), getOutputBoundaryNx()));
}

void NFFFT::update(size_t current_time_step) {
  _current_time_step = current_time_step;
  getJeJmXN();
  getJeJmXP();
  getJeJmYN();
  getJeJmYP();
}

void NFFFT::outputData() { caculateFarfield(); }

void NFFFT::getJeJmXN() {
  auto li{getOutputBoundaryStartX()};
  auto ri{getOutputBoundaryEndX()};
  auto lj{getOutputBoundaryStartY()};
  auto rj{getOutputBoundaryEndY()};
  auto nt{getCurrentTimeStep()};
  const auto& emf{getEMFInstance()};

  SpatialIndex k{0};

  for (SpatialIndex j{lj}; j < rj; ++j) {
    _jmy_xn(nt, j - lj, k) =
        -0.5 * (emf->getEz(li, j + 1, k) + emf->getEz(li, j, k));
    _jez_xn(nt, j - lj, k) =
        -0.25 * (emf->getHy(li, j, k) + emf->getHy(li - 1, j, k) +
                 emf->getHy(li, j + 1, k) + emf->getHy(lj - 1, j + 1, k));
  }
}

void NFFFT::getJeJmXP() {
  auto li{getOutputBoundaryStartX()};
  auto ri{getOutputBoundaryEndX()};
  auto lj{getOutputBoundaryStartY()};
  auto rj{getOutputBoundaryEndY()};
  auto nt{getCurrentTimeStep()};
  const auto& emf{getEMFInstance()};

  SpatialIndex k{0};

  for (SpatialIndex j{lj}; j < rj; ++j) {
    _jmy_xp(nt, j - lj, k) =
        0.5 * (emf->getEz(ri, j + 1, k) + emf->getEz(ri, j, k));
    _jez_xp(nt, j - lj, k) =
        0.25 * (emf->getHy(ri, j, k) + emf->getHy(ri - 1, j, k) +
                emf->getHy(ri, j + 1, k) + emf->getHy(ri - 1, j + 1, k));
  }
}

void NFFFT::getJeJmYN() {
  auto li{getOutputBoundaryStartX()};
  auto ri{getOutputBoundaryEndX()};
  auto lj{getOutputBoundaryStartY()};
  auto rj{getOutputBoundaryEndY()};
  auto nt{getCurrentTimeStep()};
  const auto& emf{getEMFInstance()};

  SpatialIndex k{0};

  for (SpatialIndex i{li}; i < ri; ++i) {
    _jmx_yn(nt, k, i - li) =
        0.5 * (emf->getEz(i, lj, k) + emf->getEz(i + 1, lj, k));
    _jez_yn(nt, k, i - li) =
        0.25 * (emf->getHx(i, lj, k) + emf->getHx(i, lj - 1, 0) +
                emf->getHx(i + 1, lj, k) + emf->getHx(i + 1, lj - 1, k));
  }
}

void NFFFT::getJeJmYP() {
  auto li{getOutputBoundaryStartX()};
  auto ri{getOutputBoundaryEndX()};
  auto lj{getOutputBoundaryStartY()};
  auto rj{getOutputBoundaryEndY()};
  auto nt{getCurrentTimeStep()};
  const auto& emf{getEMFInstance()};

  SpatialIndex k{0};

  for (SpatialIndex i{li}; i < ri; ++i) {
    _jmx_yp(nt, k, i - li) =
        -0.5 * (emf->getEz(i, rj, k) + emf->getEz(i + 1, rj, k));
    _jez_yp(nt, k, i - li) =
        -0.25 * (emf->getHx(i, rj, k) + emf->getHx(i, rj - 1, k) +
                 emf->getHx(i + 1, rj, k) + emf->getHx(i + 1, rj - 1, k));
  }
}

void NFFFT::caculateFarfield() {
  auto li{getOutputBoundaryStartX()};
  auto ri{getOutputBoundaryEndX()};
  auto lj{getOutputBoundaryStartY()};
  auto rj{getOutputBoundaryEndY()};
  const auto& emf{getEMFInstance()};
  auto dx{getDx()};
  auto dy{getDy()};
  auto dz{getDz()};
  auto dt{getDt()};
  auto center_x{getOutputBoundaryCenterX()};
  auto center_y{getOutputBoundaryCenterY()};
  auto center_z{getOutputBoundaryCenterZ()};
  auto phi{getFarfieldPhi()};
  auto cos_phi{cos(phi)};
  auto sin_phi{sin(phi)};
  auto dy_dz_sinth{dy * dz};
  auto dy_dz_sinphi{dy * dz * sin_phi};
  auto dy_dz_cosphi{dy * dz * cos_phi};
  auto dz_dx_sinth{dz * dx};
  auto dz_dx_sinphi{dz * dx * sin_phi};
  auto dz_dx_cosphi{dz * dx * cos_phi};

  auto num_frequency{getTotalTimeSteps() * 2 - 1};
  // TODO(franzero): refactor this
  if (_jez_xn.shape(0) < num_frequency) {
    auto size{num_frequency - _jez_xn.shape(0)};
    auto additive_shape{_jez_xn.shape()};
    additive_shape[0] = size;
    _jez_xn = xt::concatenate(
        xt::xtuple(_jez_xn, xt::zeros<double>(additive_shape)), 0);
    _jez_xp = xt::concatenate(
        xt::xtuple(_jez_xp, xt::zeros<double>(additive_shape)), 0);
    _jmy_xn = xt::concatenate(
        xt::xtuple(_jmy_xn, xt::zeros<double>(additive_shape)), 0);
    _jmy_xp = xt::concatenate(
        xt::xtuple(_jmy_xp, xt::zeros<double>(additive_shape)), 0);
    additive_shape = _jez_yn.shape();
    additive_shape[0] = size;
    _jez_yn = xt::concatenate(
        xt::xtuple(_jez_yn, xt::zeros<double>(additive_shape)), 0);
    _jez_yp = xt::concatenate(
        xt::xtuple(_jez_yp, xt::zeros<double>(additive_shape)), 0);
    _jmx_yn = xt::concatenate(
        xt::xtuple(_jmx_yn, xt::zeros<double>(additive_shape)), 0);
    _jmx_yp = xt::concatenate(
        xt::xtuple(_jmx_yp, xt::zeros<double>(additive_shape)), 0);
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

    xt::xarray<double> temp{xt::view(_jez_yn, xt::all(), k, i - li)};
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinth * exp_jk_rpr_yn;

    temp = xt::view(_jez_yp, xt::all(), k, i - li);
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinth * exp_jk_rpr_yp;

    temp = xt::view(_jmx_yn, xt::all(), k, i - li);
    l_phi +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinphi * exp_jk_rpr_yn;

    temp = xt::view(_jmx_yp, xt::all(), k, i - li);
    l_phi +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dz_dx_sinphi * exp_jk_rpr_yp;
  }

  for (SpatialIndex j{lj}; j < rj; ++j) {
    // direction: XN XP
    auto delay_y{dy * (j - center_y + 0.5) * sin_phi};
    auto exp_jk_rpr_xn{xt::exp(i1 * k_wave * (delay_y + fixed_delay_xn))};
    auto exp_jk_rpr_xp{xt::exp(i1 * k_wave * (delay_y + fixed_delay_xp))};

    xt::xarray<double> temp{xt::view(_jez_xn, xt::all(), j - lj, k)};

    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_sinth * exp_jk_rpr_xn;

    temp = xt::view(_jez_xp, xt::all(), j - lj, k);
    n_theta +=
        -xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_sinth * exp_jk_rpr_xp;

    temp = xt::view(_jmy_xn, xt::all(), j - lj, k);
    l_phi +=
        xt::fftw::fftshift(xt::fftw::fft(temp)) * dy_dz_cosphi * exp_jk_rpr_xn;

    temp = xt::view(_jmy_xp, xt::all(), j - lj, k);
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
