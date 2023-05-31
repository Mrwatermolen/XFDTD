#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "additive_source/hard_point_source.h"
#include "util/constant.h"
#include "waveform/gaussian_waveform.h"

template <typename T, typename V>
void setVectorValue(T &v, const V &value) {
  for (auto &i : v) {
    i = value;
  }
}

template <typename T>
void printP(const std::string &name, const T &t) {
  std::cout << name << ": ";
  for (auto &&e : t) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
}

int main() {
  constexpr int nc{20};
  constexpr double c{3.0e8};
  constexpr double dz{1e-3};
  constexpr double dx{1e-3};
  constexpr double dy{1e-3};
  constexpr double tau{nc * dz / (2 * c)};
  constexpr double t_0{4.5 * tau};
  constexpr size_t nz{300};
  constexpr size_t ncpml{10};
  constexpr size_t time_steps{1000};
  constexpr double sigma{std::numeric_limits<double>::epsilon() / 1000.0};

  const double dt{
      0.99 / (xfdtd::constant::C_0 *
              std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)))};
  double a{-(2 * dt / dz) / (2 * xfdtd::constant::EPSILON_0 + dt * sigma)};
  double b{-(2 * dt / dz) / (2 * xfdtd::constant::MU_0 + dt * sigma)};

  auto source{std::make_unique<xfdtd::HardPonitSource>(
      std::make_unique<xfdtd::GaussianWaveform>(1, tau, t_0),
      xfdtd::PointVector(0, 0, 150 * dz))};

  std::vector<double> ex;
  std::vector<double> cexe;
  std::vector<double> cexhy;
  std::vector<double> hy;
  std::vector<double> chyh;
  std::vector<double> chyex;
  ex.resize(nz + 1);
  setVectorValue(ex, 0);
  cexe.resize(nz + 1);
  setVectorValue(cexe, 1);
  cexhy.resize(nz + 1);
  setVectorValue(cexhy, a);
  hy.resize(nz);
  setVectorValue(hy, 0);
  chyh.resize(nz);
  setVectorValue(chyh, 1);
  chyex.resize(nz);
  setVectorValue(chyex, b);

  std::vector<double> time_array;
  time_array.resize(time_steps);
  for (size_t i{0}; i < time_steps; ++i) {
    time_array[i] = i * dt;
  }

  constexpr double e2m{xfdtd::constant::MU_0 / xfdtd::constant::EPSILON_0};
  constexpr int order{4};
  constexpr double sigma_ratio{1};
  constexpr double alpha_max{0.0};
  constexpr double alpha_min{0.01};
  constexpr double kappa_max{10};
  constexpr double sigma_max =
      sigma_ratio * (order + 1) / (150.0 * xfdtd::constant::PI * dz);

  std::array<double, ncpml> rho_e;
  std::array<double, ncpml> rho_m;
  std::array<double, ncpml> sigma_e;
  std::array<double, ncpml> sigma_m;
  std::array<double, ncpml> kappa_e;
  std::array<double, ncpml> kappa_m;
  std::array<double, ncpml> alpha_e;
  std::array<double, ncpml> alpha_m;

  std::array<double, ncpml> cpml_b_e;
  std::array<double, ncpml> cpml_a_e;
  std::array<double, ncpml> cpml_b_m;
  std::array<double, ncpml> cpml_a_m;

  std::array<double, ncpml> psi_ex;
  std::array<double, ncpml> psi_hy;

  std::array<double, ncpml> c_psi_ex;
  std::array<double, ncpml> c_psi_hy;

  for (size_t i = 0; i < ncpml; ++i) {
    rho_e[i] = (ncpml - i - 0.75) / static_cast<double>(ncpml);
    rho_m[i] = (ncpml - i - 0.25) / static_cast<double>(ncpml);
    sigma_e[i] = (sigma_max * std::pow(rho_e[i], order));
    sigma_m[i] = (sigma_max * std::pow(rho_m[i], order));
    sigma_m[i] = sigma_m[i] * e2m;
    kappa_e[i] = (1 + (kappa_max - 1) * std::pow(rho_e[i], order));
    kappa_m[i] = (1 + (kappa_max - 1) * std::pow(rho_m[i], order));
    alpha_e[i] = (alpha_min + (alpha_max - alpha_min) * (1 - rho_e[i]));
    alpha_m[i] = (alpha_min + (alpha_max - alpha_min) * (1 - rho_m[i]));
    alpha_m[i] = alpha_m[i] * e2m;

    cpml_b_e[i] = std::exp((-dt / xfdtd::constant::EPSILON_0) *
                           (sigma_e[i] / kappa_e[i] + alpha_e[i]));
    cpml_a_e[i] = (1 / dz) * (sigma_e[i] * (cpml_b_e[i] - 1)) /
                  (kappa_e[i] * (sigma_e[i] + kappa_e[i] * alpha_e[i]));

    cpml_b_m[i] = std::exp((-dt / xfdtd::constant::MU_0) *
                           (sigma_m[i] / kappa_m[i] + alpha_m[i]));
    cpml_a_m[i] = (1 / dz) * (sigma_m[i] * (cpml_b_m[i] - 1)) /
                  (kappa_m[i] * (sigma_m[i] + kappa_m[i] * alpha_m[i]));
  }

  printP("rho_e", rho_e);
  printP("sigma_e", sigma_e);
  printP("kappa_e", kappa_e);
  printP("alpha_e", alpha_e);
  printP("cpml_b_e", cpml_b_e);
  printP("cpml_a_e", cpml_a_e);

  for (size_t i = 0; i < ncpml; ++i) {
    psi_ex[i] = 0;
    psi_hy[i] = 0;
    c_psi_ex[i] = cexhy[i + 1] * dz;
    c_psi_hy[i] = chyex[i] * dz;
    cexhy[i + 1] = cexhy[i + 1] / kappa_e[i];
    chyex[i] = chyex[i] / kappa_m[i];
  }

  printP("c_psi_ex", c_psi_ex);
  printP("c_psi_hy", c_psi_hy);

  source->init(time_array);
  auto oup_dir = std::filesystem::absolute("output");
  std::cout << "Out dir: " << oup_dir.c_str() << std::endl;
  if (!std::filesystem::exists(oup_dir) &&
      !std::filesystem::is_directory(oup_dir)) {
    try {
      std::filesystem::create_directory(oup_dir);
    } catch (std::exception e) {
      std::cerr << e.what() << std::endl;
    }
  }
  for (size_t i{0}; i < time_steps; ++i) {
    // std::cout << "Time: " << i << std::endl;
    ex[150] = ex[150] + source->getValue(i) * (a * dz);
    for (size_t j{0}; j < nz; ++j) {
      hy[j] = hy[j] + chyex[j] * (ex[j + 1] - ex[j]);
    }
    for (size_t j{0}; j < ncpml; ++j) {
      psi_hy[j] = cpml_b_m[j] * psi_hy[j] + cpml_a_m[j] * (ex[j + 1] - ex[j]);
      hy[j] = hy[j] + c_psi_hy[j] * psi_hy[j];
    }
    for (size_t j{1}; j < nz; ++j) {
      ex[j] = ex[j] + cexhy[j] * (hy[j] - hy[j - 1]);
    }
    for (size_t j{0}; j < ncpml; ++j) {
      psi_ex[j] = cpml_b_e[j] * psi_ex[j] + cpml_a_e[j] * (hy[j + 1] - hy[j]);
      ex[j + 1] = ex[j + 1] + c_psi_ex[j] * psi_ex[j];
    }

    // if (i % 11 == 10) {
    //   std::stringstream ss;
    //   ss << std::setw(4) << std::setfill('0') << i << "\n";
    //   std::string s;
    //   ss >> s;
    //   std::ofstream f(std::filesystem::absolute(oup_dir) /
    //                   ("Ex-" + s + ".dat"));

    //   for (auto k{0}; k < nz; ++k) {
    //     f << k << "\t" << ex[k] << "\n";
    //   }
    //   f.close();
    // }
  }
}