#include <cstddef>
#include <memory>
#include <utility>

#include "additive_source/hard_point_source.h"
#include "boundary/perfect_match_layer.h"
#include "electromagnetic.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_2d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

void testBasic2D() {
  constexpr int nc{20};
  constexpr double c{3.0e8};
  constexpr double dx{1e-3};
  constexpr double dy{1e-3};
  constexpr double tau{nc * dx / (2 * c)};
  constexpr double t_0{4.5 * tau};
  constexpr size_t total_time_steps{1200};

  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  auto air_material{xfdtd::Material{"air", 1, 1, 0, 0, false}};
  auto free_space{xfdtd::Object{
      "free_space",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(200 * dx, 200 * dy, 0)),
      std::make_unique<xfdtd::Material>(air_material)}};

  objects.emplace_back(std::make_unique<xfdtd::Object>(std::move(free_space)));
  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "objectA",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(75 * dx, 75 * dy, 0),
                                    Eigen::Vector3d(50 * dx, 50 * dx, 0)),
      std::make_unique<xfdtd::Material>("materialA", 5.2, 1, 0, 0, false)));

  auto gaussian_waveform{xfdtd::GaussianWaveform{1, tau, t_0}};

  auto hard_point_source{xfdtd::HardPonitSource{
      std::make_unique<xfdtd::GaussianWaveform>(std::move(gaussian_waveform)),
      Eigen::Vector3d(100 * dx, 100 * dy, 0)}};
  //   sources.emplace_back(
  //       std::make_shared<xfdtd::HardPonitSource>(std::move(hard_point_source)));

  auto tfsf{xfdtd::TFSF2D{
      30, 30, xfdtd::constant::PI * 1.25, 1,
      std::make_unique<xfdtd::GaussianWaveform>(std::move(gaussian_waveform))}};

  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XN, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XP, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YN, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YP, 10));

  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(200 * dx, 200 * dy, 0)),
      xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ,
      std::filesystem::absolute("visualizing_data/2d_movie_output"), ""}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps, 10}};
  monitors.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));
  auto simulation{xfdtd::Simulation{
      dx, objects, sources, std::make_unique<xfdtd::TFSF2D>(std::move(tfsf)),
      boundaries, monitors}};
  simulation.run(total_time_steps);
  for (auto &&e : monitors) {
    e->outputData();
  }
}

int main() {
  testBasic2D();
  return 0;
}