#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <memory>
#include <xtensor/xadapt.hpp>

#include "boundary/perfect_match_layer.h"
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

void testBasic3D() {
  auto t0{std::chrono::high_resolution_clock::now()};
  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  constexpr size_t total_time_steps{2400};
  constexpr double dl{7.5e-3};

  auto domain{
      Object{"domain",
             std::make_unique<Cube>(PointVector{-30 * dl, -30 * dl, -30 * dl},
                                    PointVector{60 * dl, 60 * dl, 60 * dl}),
             Material{"air", 1, 1, 0, 0}}};
  auto scatter_object{
      Object{"a", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1),
             Material{"material_a", 3, 2, 0, 0}}};

  objects.emplace_back(std::make_shared<Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<Object>(std::move(scatter_object)));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(xfdtd::MovieMonitor{
          std::make_unique<TimeDomainFieldMonitor>(
              std::make_unique<Cube>(PointVector{-0.235, 0, -0.235},
                                     PointVector{0.45, dl, 0.45}),
              xfdtd::PlaneType::XY, xfdtd::EMComponent::EX,
              std::filesystem::path{"./visualizing_data/3d_movie_output"}, ""),
          total_time_steps, 10}));

  xt::xarray<double> frequencies{1e9};
  xt::xarray<double> theta{M_PI / 2};
  xt::xarray<double> phi{xt::linspace(0.0, 2 * M_PI, 360, false)};
  auto simulation{std::make_unique<Simulation>(
      dl, std::move(objects), std::move(sources),
      std::make_unique<TFSF3D>(
          15, 15, 15, 1, 0, 0, 0,
          std::make_unique<GaussianWaveform>(1, 2.501732435037432e-10,
                                             1.125779595766844e-09)),
      std::make_unique<NffftFd>(13, 13, 13, frequencies, theta, phi,
                                "./visualizing_data/far_feild"),
      std::move(boundaries), std::move(monitors), 0.9)};
  simulation->run(total_time_steps);
  auto dt{simulation->getDt()};
  auto waveform{
      GaussianWaveform{1, 2.501732435037432e-10, 1.125779595766844e-09}};
  waveform.init(simulation->getTimeArray());
  std::ofstream ofs{"./visualizing_data/far_feild/incident_power.dat"};
  auto incident_wave_power{
      xfdtd::dft(xt::adapt(waveform.getAllValues()), dt, frequencies)};
  for (auto i{0}; i < incident_wave_power.size(); ++i) {
    ofs << frequencies(i) << " " << std::norm(incident_wave_power(i))
        << std::endl;
  }
  ofs.close();
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << std::endl;
}

void testBandRespond() {
  auto t0{std::chrono::high_resolution_clock::now()};
  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  constexpr size_t total_time_steps{800};
  constexpr double dl{7.5e-3};

  auto domain{
      Object{"domain",
             std::make_unique<Cube>(PointVector{-30 * dl, -30 * dl, -30 * dl},
                                    PointVector{60 * dl, 60 * dl, 60 * dl}),
             Material{"air", 1, 1, 0, 0}}};
  auto scatter_object{
      Object{"a", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1),
             Material{"material_a", 3, 2, 0, 0}}};

  objects.emplace_back(std::make_shared<Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<Object>(std::move(scatter_object)));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(xfdtd::MovieMonitor{
          std::make_unique<TimeDomainFieldMonitor>(
              std::make_unique<Cube>(PointVector{-0.235, -0.235, 0},
                                     PointVector{0.45, 0.45, dl}),
              xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ,
              std::filesystem::path{"./visualizing_data/3d_movie_output"}, ""),
          total_time_steps, 10}));
  auto simulation{std::make_unique<Simulation>(
      dl, std::move(objects), std::move(sources),
      std::make_unique<TFSF3D>(
          15, 15, 15, 1, 0, 0, 0,
          std::make_unique<GaussianWaveform>(1, 2.501732435037432e-10,
                                             1.125779595766844e-09)),
      std::make_unique<xfdtd::NffftBroadBand>(13, 13, 13, M_PI, M_PI,
                                              "./visualizing_data/far_feild"),
      std::move(boundaries), std::move(monitors), 0.9)};
  simulation->run(total_time_steps);
  auto dt{simulation->getDt()};
  auto waveform{
      GaussianWaveform{1, 2.501732435037432e-10, 1.125779595766844e-09}};
  waveform.init(simulation->getTimeArray());
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << std::endl;
}

void testSphereScatter() {
  auto t0{std::chrono::high_resolution_clock::now()};
  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  constexpr size_t total_time_steps{1000};
  //   constexpr double dl{0.05};
  constexpr double dl{7.5e-3};

  auto domain{
      Object{"domain",
             std::make_unique<Cube>(PointVector{-30 * dl, -30 * dl, -30 * dl},
                                    PointVector{60 * dl, 60 * dl, 60 * dl}),
             Material{"air", 1, 1, 0, 0}}};
  auto scatter_object{
      Object{"a", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1),
             Material{"material_a", 1, 1, 1e10, 0}}};

  objects.emplace_back(std::make_shared<Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<Object>(std::move(scatter_object)));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  constexpr double dt{8.66625e-11};

  GaussianWaveform waveform{1, 2.501732435037432e-10, 1.125779595766844e-09};

  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(xfdtd::MovieMonitor{
          std::make_unique<TimeDomainFieldMonitor>(
              std::make_unique<Cube>(PointVector{-0.235, -0.235, 0},
                                     PointVector{0.45, 0.45, dl}),
              xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ,
              std::filesystem::path{"./visualizing_data/3d_movie_output"}, ""),
          total_time_steps, 10}));

  Simulation simulation{
      dl,
      objects,
      sources,
      std::make_unique<TFSF3D>(15, 15, 15, 1, 0, 0, 0,
                               std::make_unique<GaussianWaveform>(waveform)),
      std::make_unique<NffftBroadBand>(13, 13, 13, M_PI, M_PI,
                                       "./visualizing_data/far_feild"),
      boundaries,
      monitors,
      0.9};
  simulation.run(total_time_steps);
  std::vector<double> times = [&]() {
    std::vector<double> v;
    v.resize(total_time_steps * 2 - 1);
    for (int i = 0; i < v.size(); ++i) {
      v[i] = i * simulation.getDt();
    }
    return v;
  }();
  waveform.init(times);
  xt::xarray<double> tem = xt::adapt(waveform.getAllValues(),
                                     std::vector<std::size_t>{times.size()});
  std::fstream incident_wave_fft_file{"visualizing_data/incident_wave_fft.dat",
                                      std::ios::out};
  for (auto &&e : xt::abs(xt::fftw::fftshift(xt::fftw::fft(tem)))) {
    incident_wave_fft_file << e << " ";
  }
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
}

void testCoatedSphereScatter() {
  auto t0{std::chrono::high_resolution_clock::now()};
  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  constexpr size_t total_time_steps{2000};
  constexpr double dl{0.05};
  //   constexpr double dl{5.5e-3};

  auto domain{
      Object{"domain",
             std::make_unique<Cube>(PointVector{30 * dl, -30 * dl, -30 * dl},
                                    PointVector{60 * dl, 60 * dl, 60 * dl}),
             Material{"air", 1, 1, 0, 0}}};
  auto scatter_object{
      Object{"a", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1),
             Material{"material_a", 2.25, 1, 0, 0}}};
  auto coating_object{Object{
      "b", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1 + dl * 0.75),
      Material{"material_b", 4.7, 1.6, 1.6688e-3, 1.421e3}}};

  objects.emplace_back(std::make_shared<Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<Object>(std::move(coating_object)));
  objects.emplace_back(std::make_shared<Object>(std::move(scatter_object)));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  GaussianWaveform waveform{1, 0.273e-9, 0.8 * 0.273e-9};
  Simulation simulation{
      dl,
      objects,
      sources,
      std::make_unique<TFSF3D>(15, 15, 15, 1, 0, 0, 0,
                               std::make_unique<GaussianWaveform>(waveform)),
      std::make_unique<NffftBroadBand>(13, 13, 13, M_PI, M_PI,
                                       "./visualizing_data/far_feild"),
      boundaries,
      monitors,
      0.9};
  simulation.run(total_time_steps);
  std::vector<double> times = [&]() {
    std::vector<double> v;
    v.resize(total_time_steps * 2 - 1);
    for (int i = 0; i < v.size(); ++i) {
      v[i] = i * simulation.getDt();
    }
    return v;
  }();
  std::cout << simulation.getDt() << std::endl;
  waveform.init(times);
  xt::xarray<double> tem = xt::adapt(waveform.getAllValues(),
                                     std::vector<std::size_t>{times.size()});
  std::fstream incident_wave_fft_file{"visualizing_data/incident_wave_fft.dat",
                                      std::ios::out};
  for (auto &&e : xt::abs(xt::fftw::fftshift(xt::fftw::fft(tem)))) {
    incident_wave_fft_file << e << " ";
  }
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
}

void testSphereShell() {
  constexpr double dl = 1e-5;
  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  constexpr size_t total_time_steps{400};
  auto domain{
      Object{"domain",
             std::make_unique<Cube>(PointVector{-55 * dl, -55 * dl, -55 * dl},
                                    PointVector{120 * dl, 120 * dl, 160 * dl}),
             Material{"air", 1, 1, 0, 0}}};
  auto shell_object{Object{
      "shell_object", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.5e-3),
      Material{"material_a", 1.33 * 1.33, 1, 0, 0}}};

  auto shpere_objetc{Object{
      "shell_object", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.33e-3),
      Material{"material_b", 1.46 * 1.46, 2, 0, 0}}};

  objects.emplace_back(std::make_shared<Object>(domain));
  objects.emplace_back(std::make_shared<Object>(shell_object));
  objects.emplace_back(std::make_shared<Object>(shpere_objetc));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));

  xfdtd::SinusoidalWaveform sin{1, 1e12, 0};

  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(xfdtd::MovieMonitor{
          std::make_unique<TimeDomainFieldMonitor>(
              std::make_unique<Cube>(PointVector{-0.235, 0, -0.235},
                                     PointVector{0.45, dl, 0.45}),
              xfdtd::PlaneType::XY, xfdtd::EMComponent::EX,
              std::filesystem::path{"./visualizing_data/3d_movie_output"}, ""),
          total_time_steps, 10}));

  Simulation simulation{dl,
                        objects,
                        sources,
                        std::make_unique<TFSF3D>(
                            12, 12, 12, 1, 0, 0, 0,
                            std::make_unique<xfdtd::SinusoidalWaveform>(sin)),
                        nullptr,
                        boundaries,
                        monitors,
                        0.9};

  simulation.run(total_time_steps);
}

int main() { testSphereShell(); }