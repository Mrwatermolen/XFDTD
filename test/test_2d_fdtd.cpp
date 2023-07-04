#include <chrono>
#include <cstddef>
#include <fstream>
#include <memory>
#include <utility>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xmanipulation.hpp>

#include "boundary/perfect_match_layer.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "nffft/nffft.h"
#include "nffft/nffft_2d_test.h"
#include "object/object.h"
#include "shape/cylinder.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_2d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"

void testBasic2D() {
  xfdtd::SpatialIndex nx{350};
  xfdtd::SpatialIndex ny{350};
  xfdtd::SpatialIndex pml_size{10};
  xfdtd::SpatialIndex total_nx{nx + pml_size * 2};
  xfdtd::SpatialIndex total_ny{ny + pml_size * 2};
  double center_frequency{12e9};
  double max_frequency{20e9};
  double min_lambda{xfdtd::constant::C_0 / max_frequency};
  double bandwidth{2 * center_frequency};
  double dx{min_lambda / 20};
  double dy{dx};
  double tau{1.7 / (max_frequency - center_frequency)};
  double t_0{0.8 * tau};
  size_t total_time_steps{1400};
  double cylinder_x{nx * 0.5 * dx};
  double cylinder_y{ny * 0.5 * dx};
  double cylinder_radius{0.03};

  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  auto air_material{xfdtd::Material{"air", 1, 1, 0, 0, false}};
  auto pec_material{xfdtd::Material{"pec", 1, 1, 1e24, 0, false}};
  auto free_space{xfdtd::Object{
      "free_space",
      std::make_unique<xfdtd::Cube>(xfdtd::PointVector{0, 0, 0},
                                    xfdtd::PointVector{nx * dx, ny * dy, 0}),
      std::make_unique<xfdtd::Material>(air_material)}};

  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(free_space)));
  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "objectA",
      std::make_unique<xfdtd::Cylinder>(
          xfdtd::PointVector{cylinder_x, cylinder_y, 0}, cylinder_radius, 0),
      std::make_unique<xfdtd::Material>(pec_material)));

  auto cosine_modulated_gaussian_waveform{
      xfdtd::CosineModulatedGaussianWaveform{1, tau, t_0, center_frequency}};

  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XN, pml_size));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XP, pml_size));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YN, pml_size));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YP, pml_size));

  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0, 0, 0},
          xfdtd::PointVector{total_nx * dx, total_ny * dy, 0}),
      xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ,
      std::filesystem::absolute("visualizing_data/2d_movie_output"), ""}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps, 10}};
  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));
  auto simulation{xfdtd::Simulation{
      dx, objects, sources,
      std::make_unique<xfdtd::TFSF2D>(
          50, 50, xfdtd::constant::PI * 0.25, 1,
          std::make_unique<xfdtd::CosineModulatedGaussianWaveform>(
              std::move(cosine_modulated_gaussian_waveform))),
      std::make_unique<xfdtd::NFFFT2DTEST>(40, 40, 0, xfdtd::constant::PI / 2,
                                           xfdtd::constant::PI * 1.25,
                                           "visualizing_data/"),
      boundaries, monitors, 0.8}};
  auto t0{std::chrono::high_resolution_clock::now()};
  simulation.run(total_time_steps);
  for (auto &&e : monitors) {
    e->outputData();
  }

  // record incident wave fft

  std::vector<double> times = [&]() {
    std::vector<double> v;
    v.resize(total_time_steps * 2 - 1);
    for (int i = 0; i < v.size(); ++i) {
      v[i] = i * simulation.getDt();
    }
    return v;
  }();
  std::fstream incident_wave_fft_file{"visualizing_data/incident_wave_fft.dat",
                                      std::ios::out};
  auto waveform{
      xfdtd::CosineModulatedGaussianWaveform{1, tau, t_0, center_frequency}};
  waveform.init(times);
  xt::xarray<double> tem = xt::adapt(waveform.getAllValues(),
                                     std::vector<std::size_t>{times.size()});
  for (auto &&e : xt::abs(xt::fftw::fftshift(xt::fftw::fft(tem)))) {
    incident_wave_fft_file << e << " ";
  }
  incident_wave_fft_file.close();
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << std::endl;
}

int main() {
  testBasic2D();
  return 0;
}
