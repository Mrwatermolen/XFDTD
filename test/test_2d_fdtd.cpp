#include <chrono>
#include <cstddef>
#include <filesystem>
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
#include "nffft/nffft_broadband.h"
#include "object/object.h"
#include "shape/cylinder.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_2d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"

void testBasic2D() {
  constexpr xfdtd::SpatialIndex nx{350};
  constexpr xfdtd::SpatialIndex ny{300};
  constexpr xfdtd::SpatialIndex pml_size{10};
  constexpr xfdtd::SpatialIndex total_nx{nx + pml_size * 2};
  constexpr xfdtd::SpatialIndex total_ny{ny + pml_size * 2};
  constexpr double center_frequency{12e9};
  constexpr double max_frequency{20e9};
  constexpr double min_lambda{3e8 / max_frequency};
  constexpr double bandwidth{2 * center_frequency};
  constexpr double dx{min_lambda / 20};
  constexpr double dy{dx};
  constexpr double tau{1.7 / (max_frequency - center_frequency)};
  constexpr double t_0{0.8 * tau};
  constexpr size_t total_time_steps{1400};
  constexpr double cylinder_x{nx * 0.5 * dx};
  constexpr double cylinder_y{ny * 0.5 * dx};
  constexpr double cylinder_radius{0.03};
  const std::filesystem::path output_dir{
      "./visualizing_data/data/2d_cylinder_scatter/"};

  auto objects{xfdtd::ObjectArray{}};
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
          xfdtd::Axis::Z, xfdtd::PointVector{cylinder_x, cylinder_y, 0},
          cylinder_radius, 0),
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

  const std::filesystem::path movie_dir{output_dir / "movie_ez"};
  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<xfdtd::Cube>(xfdtd::PointVector{0, 0, 0},
                                    xfdtd::PointVector{nx * dx, ny * dy, 0}),
      xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ, movie_dir, ""}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps, 100}};

  auto simulation{xfdtd::Simulation{
      dx, objects, boundaries,
      std::make_unique<xfdtd::TFSF2D>(
          50, 50, xfdtd::constant::PI * 0.25, 1,
          std::make_unique<xfdtd::CosineModulatedGaussianWaveform>(
              std::move(cosine_modulated_gaussian_waveform))),
      std::make_unique<xfdtd::NffftBroadBand>(
          40, 40, 0, xfdtd::constant::PI / 2, xfdtd::constant::PI * 1.25,
          output_dir / "nffft"),
      0.8}};
  simulation.addMonitor(
      std::make_unique<xfdtd::MovieMonitor>(std::move(movie_monitor)));

  auto t0{std::chrono::high_resolution_clock::now()};
  simulation.run(total_time_steps);
  simulation.outputTFSFIncidentWaveFastFourierTransform(
      output_dir / "tfsf_incident_wave_fft");
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << '\n';
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << '\n';
}

int main() {
  testBasic2D();
  return 0;
}
