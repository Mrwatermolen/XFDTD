#include <chrono>
#include <filesystem>
#include <memory>
#include <utility>

#include "boundary/perfect_match_layer.h"
#include "electromagnetic_field/electromagnetic_field.h"
#include "helper.h"
#include "material/material.h"
#include "monitor/field_monitor.h"
#include "monitor/movie_monitor.h"
#include "nffft/nffft_fd.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/sphere.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_3d.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

using xfdtd::Boundary;
using xfdtd::Cube;
using xfdtd::GaussianWaveform;
using xfdtd::Material;
using xfdtd::Monitor;
using xfdtd::NffftFd;
using xfdtd::Object;
using xfdtd::PML;
using xfdtd::PointVector;
using xfdtd::Simulation;
using xfdtd::Sphere;
using xfdtd::TFSF3D;

void pecSphereBistaticRCS() {
  // Grid
  constexpr double dl = 7.5e-3;

  auto objects{std::vector<std::shared_ptr<Object>>{}};
  auto domain_origin_point{-40 * dl};
  auto domain_size{80 * dl};
  objects.emplace_back(std::make_shared<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{domain_origin_point, domain_origin_point,
                      domain_origin_point},
          PointVector{domain_size, domain_size, domain_size}),
      Material{"air", 1, 1, 0, 0}));

  constexpr double radius{0.1};
  objects.emplace_back(std::make_shared<Object>(
      "scatter", std::make_unique<Sphere>(PointVector{0, 0, 0}, radius),
      Material{"pec", 1, 1, 1e10, 0}));

  auto boundaries{std::vector<std::shared_ptr<Boundary>>{}};
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
      std::filesystem::path{
          "./visualizing_data/data/pec_sphere_bistatic_rcs"})};

  auto monitor{std::make_unique<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::TimeDomainFieldMonitor>(
          std::make_unique<Cube>(
              PointVector{domain_origin_point, 0, domain_origin_point},
              PointVector{domain_size, 0, domain_size}),
          xfdtd::PlaneType::ZX, xfdtd::EMComponent::EZ,
          "./visualizing_data/data/pec_sphere_bistatic_rcs/movie", ""),
      500, 10)};

  constexpr size_t total_time_steps{500};
  Simulation simulation(dl, objects, boundaries, std::move(tfsf),
                        std::move(bistatic_rcs), 1);
  // simulation.addMonitor(std::move(monitor));
  simulation.run(total_time_steps);
}

int main() {
  auto duration{xfdtd_example::timeSomething(pecSphereBistaticRCS)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::milliseconds>(duration)};
  std::cout << "It costs " << duration_in_seconds.count() << " seconds or "
            << duration_in_milliseconds.count()
            << " milliseconds to run PEC Sphere Bistatic RCS"
            << "\n";
}