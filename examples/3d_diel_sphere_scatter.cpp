#include <filesystem>
#include <memory>
#include <xtensor/xnpy.hpp>

#include "boundary/perfect_match_layer.h"
#include "helper.h"
#include "material/material.h"
#include "nffft/nffft_fd.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/sphere.h"
#include "simulation/simulation.h"
#include "tfsf/tfsf_3d.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

void dielectricSphereScatter() {
  constexpr double dl{7.5e-3};

  auto simulation{xfdtd::Simulation{dl}};

  auto sphere{xfdtd::Material{"sphere", 3, 2, 0, 0}};

  constexpr double radius{0.1};
  constexpr double domain_size{2 * radius + 18 * dl};
  constexpr int num_cells{static_cast<int>(domain_size / dl)};
  auto domain_shape{xfdtd::Cube{
      xfdtd::PointVector{-domain_size / 2, -domain_size / 2, -domain_size / 2},
      xfdtd::PointVector{domain_size, domain_size, domain_size}}};
  auto sphere_shape{xfdtd::Sphere{xfdtd::PointVector{0, 0, 0}, radius}};

  auto domain{xfdtd::Object{
      "domain", std::make_unique<xfdtd::Cube>(std::move(domain_shape)),
      xfdtd::Material::createAir()}};
  auto sphere_object{xfdtd::Object{
      "sphere", std::make_unique<xfdtd::Sphere>(std::move(sphere_shape)),
      std::move(sphere)}};

  simulation.addObject(std::make_shared<xfdtd::Object>(std::move(domain)));
  simulation.addObject(
      std::make_shared<xfdtd::Object>(std::move(sphere_object)));

  constexpr size_t pml_thickness{8};
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XN, pml_thickness));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XP, pml_thickness));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YN, pml_thickness));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YP, pml_thickness));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZN, pml_thickness));
  simulation.addBoundary(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::ZP, pml_thickness));

  constexpr size_t n{20};
  constexpr double lambda_min{dl * n};
  constexpr double f_max{3e8 / lambda_min};
  constexpr double tau{lambda_min / 6e8};
  constexpr double t_0{4.5 * tau};
  auto waveform{xfdtd::GaussianWaveform{tau, t_0}};

  constexpr size_t inject_index{pml_thickness + 10};
  constexpr double incident_theta{0};
  constexpr double incident_phi{0};
  constexpr double polarization{0};
  constexpr double incident_amplitude{1};
  auto tfsf_ptr{std::make_shared<xfdtd::TFSF3D>(
      inject_index, inject_index, inject_index, incident_amplitude,
      incident_theta, incident_theta, polarization,
      std::make_unique<xfdtd::GaussianWaveform>(waveform))};

  simulation.addTFSFSource(tfsf_ptr);

  auto frequencies{1e9};
  auto phi{xt::xarray<double>{0}};
  auto theta{xt::linspace<double>(0, xfdtd::constant::PI, 181)};
  auto nffft_ptr{std::make_shared<xfdtd::NffftFd>(
      inject_index - 3, inject_index - 3, inject_index - 3, frequencies, theta,
      phi,
      std::filesystem::path{
          "./visualizing_data/data/3d_dielectric_sphere_scatter"})};

  simulation.addNFFFT(nffft_ptr);

  simulation.run(2000);

  auto incident_fft{
      tfsf_ptr->getIncidentWaveFourierTransform(xt::xarray<double>{1e9})};

  xt::dump_npy(
      "./visualizing_data/data/3d_dielectric_sphere_scatter/incident_fft.npy",
      incident_fft);
}

int main() {
  auto duration{xfdtd_example::timeSomething(dielectricSphereScatter)};
  std::cout
      << "It costs " << duration.count() << " seconds or "
      << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
      << " milliseconds to run Dielectric Sphere Scatter."
      << "\n";
}