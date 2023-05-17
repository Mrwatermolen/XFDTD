#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "electromagnetic.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"
#include "waveform/sinusoidal_waveform.h"

void testBasic() {
  constexpr int nc{20};
  constexpr double c{3.0e8};
  constexpr double dz{1e-3};
  constexpr double tau{nc * dz / (2 * c)};
  constexpr double t_0{4.5 * tau};
  constexpr size_t total_time_steps{1200};

  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  auto air_material{xfdtd::Material{"vaccum", 1, 1, 0, 0, false}};
  auto free_space{xfdtd::Object{
      "free_space",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(0, 0, 400 * dz)),
      std::make_unique<xfdtd::Material>(air_material)}};
  objects.emplace_back(std::make_shared<xfdtd::Object>(free_space));
  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "objectA",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 100 * dz),
                                    Eigen::Vector3d(0, 0, 150 * dz)),
      std::make_unique<xfdtd::Material>("materialA", 5, 1.4, 0, 0, false)));

  auto gaussian_waveform{xfdtd::GaussianWaveform{1, tau, t_0}};
  auto gaussian_point_source{xfdtd::HardPonitSource{
      std::make_unique<xfdtd::GaussianWaveform>(std::move(gaussian_waveform)),
      Eigen::Vector3d(0, 0, 250 * dz)}};

  sources.emplace_back(std::make_shared<xfdtd::HardPonitSource>(
      std::move(gaussian_point_source)));

  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZN, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZP, 10));

  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, -10 * dz),
                                    Eigen::Vector3d(0, 0, 410 * dz)),
      xfdtd::PlaneType::ZX, xfdtd::EMComponent::EX,
      std::filesystem::absolute("movie_monitor"), "ex.dat"}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps}};
  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));
  //   monitors.push_back(
  //       std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));
  //   xfdtd::MovieMonitor movie_monitor{
  //       std::move(std::make_shared<xfdtd::TimeDomainFieldMonitor>(
  //           std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
  //                                         Eigen::Vector3d(0, 0, 400 * dz)),
  //           xfdtd::PlaneType::ZX, xfdtd::EMComponent::EX,
  //           std::filesystem::absolute("movie_monitor"), "")),
  //       total_time_steps};

  //   monitors.emplace_back(std::move(movie_monitor));
  //   monitors.emplace_back(std::make_shared<xfdtd::TimeDomainFieldMonitor>(
  //       std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
  //                                     Eigen::Vector3d(0, 0, 400 * dz)),
  //       xfdtd::PlaneType::ZX, xfdtd::EMComponent::EX,
  //       std::filesystem::absolute("monitor"), "ex.dat"));

  auto simulation{
      xfdtd::Simulation(dz, objects, sources, boundaries, monitors)};
  simulation.run(total_time_steps);
  for (auto& m : monitors) {
    m->outputData();
  }
}

int main() { testBasic(); }
