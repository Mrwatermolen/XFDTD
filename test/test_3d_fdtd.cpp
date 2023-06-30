#include <cmath>
#include <filesystem>
#include <memory>

#include "boundary/perfect_match_layer.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/sphere.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_3d.h"
#include "waveform/gaussian_waveform.h"

void testBasic3D() {
  using xfdtd::Cube;
  using xfdtd::GaussianWaveform;
  using xfdtd::Material;
  using xfdtd::Object;
  using xfdtd::PML;
  using xfdtd::PointVector;
  using xfdtd::Simulation;
  using xfdtd::Sphere;
  using xfdtd::TFSF3D;
  using xfdtd::TimeDomainFieldMonitor;

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
             Material{"pec", 1, 1, 1e20, 0}}};

  objects.emplace_back(std::make_shared<Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<Object>(std::move(scatter_object)));

  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XN, 10));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::XP, 10));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YN, 10));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::YP, 10));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZN, 10));
  boundaries.emplace_back(std::make_shared<PML>(xfdtd::Orientation::ZP, 10));

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
          15, 15, 15, 1, M_PI / 2, M_PI / 4, M_PI / 4,
          std::make_unique<GaussianWaveform>(1, 2.501732435037432e-10,
                                             1.125779595766844e-09)),
      nullptr, std::move(boundaries), std::move(monitors), 0.9)};
  simulation->run(total_time_steps);
  auto t1{std::chrono::high_resolution_clock::now()};
  std::cout
      << "Simulation takes "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " msec" << std::endl;
  std::cout << "Simulation takes "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s" << std::endl;
}

int main() { testBasic3D(); }