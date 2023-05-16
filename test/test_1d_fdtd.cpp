#include <iostream>
#include <memory>
#include <vector>

#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "object/object.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"
#include "waveform/sinusoidal_waveform.h"

void test() {
  constexpr int nc{20};
  constexpr double c{3.0e8};
  constexpr double dz{1e-3};
  constexpr double tau{nc * dz / (2 * c)};
  constexpr double t_0{4.5 * tau};

  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};

  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "free_space",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(0, 0, 400 * dz)),
      std::make_unique<xfdtd::Material>("vaccum", 1, 1, 0, 0, false)));

  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "objectA",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 100 * dz),
                                    Eigen::Vector3d(0, 0, 150 * dz)),
      std::make_unique<xfdtd::Material>("materialA", 5, 1.4, 0, 0, false)));

  sources.emplace_back(std::make_shared<xfdtd::HardPonitSource>(
      std::make_unique<xfdtd::GaussianWaveform>(1, tau, t_0),
      Eigen::Vector3d(0, 0, 250 * dz)));

  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZN, 10));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZP, 10));

  auto simulation{xfdtd::Simulation(dz, objects, sources, boundaries)};
  simulation.run(1200);
}

int main() { test(); }
