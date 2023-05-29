#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_1d.h"
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
      std::make_unique<xfdtd::Cube>(xfdtd::PointVector(0, 0, 0),
                                    xfdtd::PointVector(0, 0, 300 * dz)),
      std::make_unique<xfdtd::Material>(air_material)}};
  objects.emplace_back(std::make_shared<xfdtd::Object>(free_space));
  //   objects.emplace_back(std::make_shared<xfdtd::Object>(
  //       "objectA",
  //       std::make_unique<xfdtd::Cube>(xfdtd::PointVector(0, 0, 100 * dz),
  //                                     xfdtd::PointVector(0, 0, 150 * dz)),
  //       std::make_unique<xfdtd::Material>("materialA", 2.2, 1, 0, 0,
  //       false)));

  auto gaussian_waveform{xfdtd::GaussianWaveform{1, tau, t_0}};

  auto gaussian_point_source{xfdtd::HardPonitSource{
      std::make_unique<xfdtd::GaussianWaveform>(std::move(gaussian_waveform)),
      xfdtd::PointVector(0, 0, 50 * dz)}};
  sources.emplace_back(std::make_shared<xfdtd::HardPonitSource>(
      std::move(gaussian_point_source)));

  //   auto tfsf{xfdtd::TFSF1D{
  //       30, xfdtd::constant::PI, -1,
  //       std::make_unique<xfdtd::GaussianWaveform>(std::move(gaussian_waveform))}};

  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZN, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZP, 10));

  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<xfdtd::Cube>(xfdtd::PointVector(0, 0, -10 * dz),
                                    xfdtd::PointVector(0, 0, 310 * dz)),
      xfdtd::PlaneType::ZX, xfdtd::EMComponent::EX,
      std::filesystem::absolute("visualizing_data/1d_movie_monitor"), ""}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps, 20}};
  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));

  auto simulation{
      xfdtd::Simulation(dz, objects, sources, boundaries, monitors)};
  simulation.run(total_time_steps);
  for (auto& m : monitors) {
    m->outputData();
  }
}

int main() { testBasic(); }
