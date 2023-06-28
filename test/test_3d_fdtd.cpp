#include <memory>

#include "additive_source/hard_point_source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "shape/sphere.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_3d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

using xfdtd::Cube;
using xfdtd::GaussianWaveform;
using xfdtd::HardPonitSource;
using xfdtd::Material;
using xfdtd::Object;
using xfdtd::PML;
using xfdtd::PointVector;
using xfdtd::Sphere;

void testBasic3D() {
  constexpr double dx{0.0075};
  constexpr double dy{0.0075};
  constexpr double dz{0.0075};
  constexpr float cfl{0.99};
  constexpr size_t total_time_steps{1000};

  xfdtd::ObjectArray object_array;
  xfdtd::SourceArray source_array;
  xfdtd::BoundaryArray boundary_array;
  xfdtd::MonitorArray monitor_array;

  auto air{Material{"air", 1, 1, 0, 0}};
  auto domain{Cube{PointVector{-1, -1, -1}, PointVector{2, 2, 2}}};
  auto free_space{
      Object{"free_space",
             std::make_unique<Cube>(PointVector{-0.1750, -0.1750, -0.1750},
                                    PointVector{0.3500, 0.3500, 0.3500}),
             air}};

  auto pec{Material{"pec", 1, 1, 1e24, 0}};
  auto sphere{Sphere{PointVector{0.4, 0.4, 0}, 0.2}};
  auto object_a{
      Object{"a", std::make_unique<Sphere>(PointVector{0, 0, 0}, 0.1), pec}};

  object_array.emplace_back(std::make_shared<Object>(free_space));
  object_array.emplace_back(std::make_shared<Object>(object_a));

  double tau{30 * cfl * 0.5 / xfdtd::constant::C_0};
  auto gaussian_waveform{GaussianWaveform{1, tau, 4.5 * tau}};

  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 8));
  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 8));
  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 8));
  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 8));
  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 8));
  boundary_array.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 8));

  auto monitor{xfdtd::TimeDomainFieldMonitor{
      std::make_unique<Cube>(PointVector{-1.5, -1.5, 0}, PointVector{3, 3, dz}),
      xfdtd::PlaneType::XY, xfdtd::EMComponent::EZ,
      std::filesystem::absolute("visualizing_data/3d_movie_output"), ""}};
  auto movie_monitor{xfdtd::MovieMonitor{
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(std::move(monitor)),
      total_time_steps, 10}};
  monitor_array.emplace_back(
      std::make_shared<xfdtd::MovieMonitor>(std::move(movie_monitor)));
  auto t0{std::chrono::high_resolution_clock::now()};
  auto simulation{xfdtd::Simulation{
      dx, object_array, source_array,
      std::make_unique<xfdtd::TFSF3D>(
          15, 15, 15, 1000000, xfdtd::constant::PI / 2,
          xfdtd::constant::PI * 0.25, xfdtd::constant::PI * 0.25,
          std::make_unique<xfdtd::GaussianWaveform>(gaussian_waveform)),
      nullptr, boundary_array, monitor_array, cfl}};
  simulation.checkRun(total_time_steps);
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << std::endl;
  int a;
  std::cin >> a;
  simulation.checkRun(1);
}

int main() { testBasic3D(); }