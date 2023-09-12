#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <utility>
#include <xtensor/xadapt.hpp>

#include "boundary/perfect_match_layer.h"
#include "material/dispersive_material.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "nffft/nffft_broadband.h"
#include "nffft/nffft_fd.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/sphere.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_3d.h"
#include "util/dft.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"
#include "waveform/gaussian_waveform.h"
#include "waveform/sinusoidal_waveform.h"

using xfdtd::Cube;
using xfdtd::GaussianWaveform;
using xfdtd::Material;
using xfdtd::NffftBroadBand;
using xfdtd::NffftFd;
using xfdtd::Object;
using xfdtd::PML;
using xfdtd::PointVector;
using xfdtd::Simulation;
using xfdtd::Sphere;
using xfdtd::TFSF3D;
using xfdtd::TimeDomainFieldMonitor;

void testPECSphereMonostaticRCS() {
  // Grid
  constexpr double dl = 2e-2;

  auto objects{xfdtd::ObjectArray{}};
  auto domain_orign_point{-40 * dl};
  auto domain_size{80 * dl};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_orign_point, domain_orign_point,
                      domain_orign_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));

  constexpr double radius{0.5};
  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      Material{"pec", 1, 1, 1e10, 0}));

  auto boundaries{xfdtd::BoundaryArray{}};
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  constexpr int n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  constexpr size_t tfsf_boundary_index{15};
  auto tfsf{std::make_unique<TFSF3D>(
      tfsf_boundary_index, tfsf_boundary_index, tfsf_boundary_index,
      incident_amplitude, incident_theta, incident_phi, polarization,
      std::make_unique<GaussianWaveform>(1, tau, t_0))};

  constexpr size_t output_boundary_index{12};
  auto mono_rcs{std::make_unique<NffftBroadBand>(
      output_boundary_index, output_boundary_index, output_boundary_index,
      incident_theta + M_PI, incident_phi + M_PI,
      std::filesystem::path{"./visualizing_data/mie"})};

  constexpr size_t total_time_steps{1200};
  Simulation simulation(dl, objects, boundaries, std::move(tfsf),
                        std::move(mono_rcs),  0.98);
  simulation.run(total_time_steps);
}

void testPECSphereBistaticRCS() {
  // Grid
  constexpr double dl = 7.5e-3;

  auto objects{xfdtd::ObjectArray{}};
  auto domain_orign_point{-40 * dl};
  auto domain_size{80 * dl};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_orign_point, domain_orign_point,
                      domain_orign_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));

  constexpr double radius{0.1};
  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      Material{"pec", 1, 1, 1e10, 0}));

  auto boundaries{xfdtd::BoundaryArray{}};
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 12));

  constexpr int n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  constexpr size_t tfsf_boundary_index{27};
  auto tfsf{std::make_unique<TFSF3D>(
      tfsf_boundary_index, tfsf_boundary_index, tfsf_boundary_index,
      incident_amplitude, incident_theta, incident_phi, polarization,
      std::make_unique<GaussianWaveform>(1, tau, t_0))};

  constexpr size_t output_boundary_index{17};
  auto frequencies{1e9};
  auto phi{xt::xarray<double>{0}};
  auto theta{xt::linspace(0.0, M_PI, 180, false)};
  auto bistatic_rcs{std::make_unique<NffftFd>(
      output_boundary_index, output_boundary_index, output_boundary_index,
      frequencies, theta, phi,
      std::filesystem::path{"./visualizing_data/bistatic_rcs"})};

  constexpr size_t total_time_steps{1000};
  Simulation simulation(dl, objects, boundaries, std::move(tfsf),
                        std::move(bistatic_rcs), 0.98);
  simulation.run(total_time_steps);
  auto waveform{xfdtd::CosineModulatedGaussianWaveform{1, tau, t_0, 1e9}};
  auto dt{simulation.getDt()};
  waveform.init(simulation.getTimeArray());
  std::ofstream ofs{"./visualizing_data/bistatic_rcs/incident_power.dat"};
  std::cout << simulation.getDt();
  auto incident_wave_power{
      xfdtd::dft(xt::adapt(waveform.getAllValues()), dt, frequencies)};
  for (auto i{0}; i < incident_wave_power.size(); ++i) {
    ofs << std::norm(incident_wave_power(i)) << "\t";
  }
  ofs.close();
}

void testLorentzSphereMonostaticRCS() {
  constexpr double radius{15e-9};
  constexpr double dl{radius / 30};
  constexpr double domain_size{radius * 2.5};
  constexpr double domain_orign_point{-domain_size / 2};
  constexpr double domain_cell_x{domain_size / dl};
  auto objects{xfdtd::ObjectArray{}};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_orign_point, domain_orign_point,
                      domain_orign_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));
  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      std::make_unique<xfdtd::LorentzMedium>("random", 2.25, 1, 4e16,
                                             0.28e16)));
  auto boundaries{xfdtd::BoundaryArray{}};
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 12));

  constexpr int n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  constexpr size_t tfsf_boundary_index{22};
  auto tfsf{std::make_unique<TFSF3D>(
      tfsf_boundary_index, tfsf_boundary_index, tfsf_boundary_index,
      incident_amplitude, incident_theta, incident_phi, polarization,
      std::make_unique<GaussianWaveform>(1, tau, t_0))};

  constexpr size_t output_boundary_index{15};
  auto mono_rcs{std::make_unique<NffftBroadBand>(
      output_boundary_index, output_boundary_index, output_boundary_index,
      incident_theta + M_PI, incident_phi + M_PI,
      std::filesystem::path{"./visualizing_data/mono_rcs"})};

  Simulation simulation(dl, objects, boundaries, std::move(tfsf),
                        std::move(mono_rcs));
  constexpr size_t total_time_steps{1000};
  simulation.run(total_time_steps);
}

void testDrudeSphereMonostaticRCS() {
  constexpr double radius{3.75e-3};
  constexpr double dl{radius / 30};
  constexpr double domain_size{radius * 2.5};
  constexpr double domain_orign_point{-domain_size / 2};
  constexpr double domain_cell_x{domain_size / dl};
  auto objects{xfdtd::ObjectArray{}};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_orign_point, domain_orign_point,
                      domain_orign_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));

  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      std::make_unique<xfdtd::DrudeMedium>("random", 1, 1.8e11, 2e10)));
  auto boundaries{xfdtd::BoundaryArray{}};
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 12));

  constexpr int n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  constexpr size_t tfsf_boundary_index{22};
  auto tfsf{std::make_unique<TFSF3D>(
      tfsf_boundary_index, tfsf_boundary_index, tfsf_boundary_index,
      incident_amplitude, incident_theta, incident_phi, polarization,
      std::make_unique<GaussianWaveform>(1, tau, t_0))};

  constexpr size_t output_boundary_index{15};
  auto mono_rcs{std::make_unique<NffftBroadBand>(
      output_boundary_index, output_boundary_index, output_boundary_index,
      incident_theta + M_PI, incident_phi + M_PI,
      std::filesystem::path{"./visualizing_data/mono_rcs"})};

  Simulation simulation(dl, objects, boundaries, std::move(tfsf),
                        std::move(mono_rcs));
  constexpr size_t total_time_steps{1000};
  simulation.run(total_time_steps);
}

void testDebySphereMonostaticRCS() {
  constexpr double radius{0.25};
  constexpr double dl{radius / 30};
  constexpr double domain_size{radius * 2.5};
  constexpr double domain_orign_point{-domain_size / 2};
  constexpr double domain_cell_x{domain_size / dl};
  auto objects{xfdtd::ObjectArray{}};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_orign_point, domain_orign_point,
                      domain_orign_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));

  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      std::make_unique<xfdtd::DebyMedium>("random", 1.16, 1.01, 6.497e-10,
                                          2.95e-4)));

  auto boundaries{xfdtd::BoundaryArray{}};
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 12));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 12));

  constexpr int n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  constexpr size_t tfsf_boundary_index{22};
  auto tfsf{std::make_unique<TFSF3D>(
      tfsf_boundary_index, tfsf_boundary_index, tfsf_boundary_index,
      incident_amplitude, incident_theta, incident_phi, polarization,
      std::make_unique<GaussianWaveform>(1, tau, t_0))};

  constexpr size_t output_boundary_index{15};
  auto mono_rcs{std::make_unique<NffftBroadBand>(
      output_boundary_index, output_boundary_index, output_boundary_index,
      incident_theta + M_PI, incident_phi + M_PI,
      std::filesystem::path{"./visualizing_data/mono_rcs"})};

  Simulation simulation{dl, objects, boundaries, std::move(tfsf),
                        std::move(mono_rcs)};
  constexpr size_t total_time_steps{1000};
  simulation.run(total_time_steps);
}

int main() {
  auto t0{std::chrono::high_resolution_clock::now()};
  testLorentzSphereMonostaticRCS();
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec "
      << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
      << " s" << std::endl;
}