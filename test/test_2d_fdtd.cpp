#include <cstddef>
#include <memory>
#include <utility>

#include "boundary/perfect_match_layer.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "object/object.h"
#include "shape/cylinder.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_2d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"

void testBasic2D() {
  xfdtd::SpatialIndex nx{150};
  xfdtd::SpatialIndex ny{150};
  xfdtd::SpatialIndex pml_size{10};
  xfdtd::SpatialIndex total_nx{nx + pml_size * 2};
  xfdtd::SpatialIndex total_ny{ny + pml_size * 2};
  double center_frequency{4e9};
  double max_frequency{8e9};
  double min_lambda{xfdtd::constant::C_0 / max_frequency};
  double bandwidth{2 * center_frequency};
  double dx{min_lambda / 20};
  double dy{dx};
  double tau{1.5 / (max_frequency - center_frequency)};
  double t_0{0.8 * tau};
  size_t total_time_steps{640};
  double cylinder_x{nx * 0.5 * dx};
  double cylinder_y{ny * 0.5 * dx};
  double cylinder_radius{0.03};

  auto objects{xfdtd::ObjectArray{}};
  auto sources{xfdtd::SourceArray{}};
  auto boundaries{xfdtd::BoundaryArray{}};
  auto monitors{xfdtd::MonitorArray{}};

  auto air_material{xfdtd::Material{"air", 1, 1, 0, 0, false}};
  auto pec_material{xfdtd::Material{"pec", 2, 1, 1e24, 0, false}};
  auto free_space{xfdtd::Object{
      "free_space",
      std::make_unique<xfdtd::Cube>(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(nx * dx, ny * dy, 0)),
      std::make_unique<xfdtd::Material>(air_material)}};

  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(free_space)));
  objects.emplace_back(std::make_shared<xfdtd::Object>(
      "objectA",
      std::make_unique<xfdtd::Cylinder>(
          xfdtd::PointVector(cylinder_x, cylinder_y, 0), cylinder_radius, 0),
      std::make_unique<xfdtd::Material>(pec_material)));

  auto cosine_modulated_gaussian_waveform{
      xfdtd::CosineModulatedGaussianWaveform{1, tau, t_0, center_frequency}};

  auto tfsf{
      xfdtd::TFSF2D{20, 20, xfdtd::constant::PI * 1.25, 1,
                    std::make_unique<xfdtd::CosineModulatedGaussianWaveform>(
                        std::move(cosine_modulated_gaussian_waveform))}};

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
          Eigen::Vector3d(0, 0, 0),
          Eigen::Vector3d(total_nx * dx, total_ny * dy, 0)),
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
