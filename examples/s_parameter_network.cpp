#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "boundary/perfect_match_layer.h"
#include "helper.h"
#include "lumped_element/resistor.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "monitor/voltage_monitor.h"
#include "network/network.h"
#include "network/port.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

void quarterWaveTransformer() {
  // define material
  auto air{xfdtd::Material{"air", 1, 1, 0, 0}};
  auto pec{xfdtd::Material{"pec", 1, 1, 1e10, 0}};

  // define grid space
  constexpr double size{2e-4};

  // define object

  // domain
  auto domain_origin_point{xfdtd::PointVector{-5 * size, -5 * size, 0 * size}};
  auto domain_size{xfdtd::PointVector{60 * size, 100 * size, 15 * size}};
  auto domain{xfdtd::Object{
      "domain", std::make_unique<xfdtd::Cube>(domain_origin_point, domain_size),
      std::make_unique<xfdtd::Material>(air)}};

  auto substrate_origin_point{xfdtd::PointVector{0, 0, 0}};
  constexpr double x{50 * size};
  constexpr double y{90 * size};
  constexpr double z{5 * size};
  auto substrate_size{xfdtd::PointVector{x, y, z}};
  auto substrate_shape{xfdtd::Cube{substrate_origin_point, substrate_size}};
  auto substrate_material{xfdtd::Material{"substrate", 4.6, 1, 0, 0}};
  auto substrate{
      xfdtd::Object{"substrate", std::make_unique<xfdtd::Cube>(substrate_shape),
                    std::make_unique<xfdtd::Material>(substrate_material)}};

  // 50 ohm microstrip line
  constexpr double msl_50_width{9 * size};
  constexpr double msl_50_length{20 * size};
  auto msl_50_origin_point{xfdtd::PointVector{4e-3, 0, z}};
  auto msl_50_size{xfdtd::PointVector{msl_50_width, msl_50_length, 0}};
  auto msl_50{xfdtd::ObjectPlane{
      "msl_50", std::make_unique<xfdtd::Cube>(msl_50_origin_point, msl_50_size),
      std::make_unique<xfdtd::Material>(pec)}};

  // 70.7 ohm microstrip line
  constexpr double msl_70_7_width{5 * size};
  constexpr double msl_70_7_length{50 * size};
  auto msl_70_7_origin_point{xfdtd::PointVector{4.4e-3, msl_50_length, z}};
  auto msl_70_7_size{xfdtd::PointVector{msl_70_7_width, msl_70_7_length, 0}};
  auto msl_70_7{xfdtd::ObjectPlane{
      "msl_70_7",
      std::make_unique<xfdtd::Cube>(msl_70_7_origin_point, msl_70_7_size),
      std::make_unique<xfdtd::Material>(pec)}};

  // 100 ohm microstrip line
  constexpr double msl_100_width{2 * size};
  constexpr double msl_100_length{20 * size};
  auto msl_100_origin_point{
      xfdtd::PointVector{4.8e-3, msl_50_length + msl_70_7_length, z}};
  auto msl_100_size{xfdtd::PointVector{msl_100_width, msl_100_length, 0}};
  auto msl_100{xfdtd::ObjectPlane{
      "msl_100",
      std::make_unique<xfdtd::Cube>(msl_100_origin_point, msl_100_size),
      std::make_unique<xfdtd::Material>(pec)}};

  auto objects{std::vector<std::shared_ptr<Object>>{}};
  objects.emplace_back(std::make_unique<xfdtd::Object>(domain));
  objects.emplace_back(std::make_unique<xfdtd::Object>(substrate));
  objects.emplace_back(std::make_unique<xfdtd::ObjectPlane>(msl_50));
  objects.emplace_back(std::make_unique<xfdtd::ObjectPlane>(msl_70_7));
  objects.emplace_back(std::make_unique<xfdtd::ObjectPlane>(msl_100));

  // define lumped element
  auto lumped_elements{std::vector<std::shared_ptr<xfdtd::LumpedElement>>()};
  constexpr double l_min{size * 20};
  auto tau{l_min / 6e8};
  auto t_0{4.5 * tau};
  auto waveform{xfdtd::GaussianWaveform{1, tau, t_0}};
  auto source{xfdtd::VoltageSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4e-3, 0, 0},
          xfdtd::PointVector{msl_50_width, 2 * size, z}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::GaussianWaveform>(std::move(waveform))}};
  auto resistor{
      xfdtd::Resistor{std::make_unique<xfdtd::Cube>(
                          xfdtd::PointVector{4.8e-3, y - 2 * size, 0},
                          xfdtd::PointVector{msl_100_width, 2 * size, z}),
                      xfdtd::Axis::Z, 100}};
  lumped_elements.emplace_back(std::make_shared<xfdtd::VoltageSource>(source));
  lumped_elements.emplace_back(std::make_shared<xfdtd::Resistor>(resistor));

  // define boundary
  auto boundaries{std::vector<std::shared_ptr<Boundary>>{}};
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XN, 8));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::XP, 8));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YN, 8));
  boundaries.emplace_back(
      std::make_shared<xfdtd::PML>(xfdtd::Orientation::YP, 8));

  // define monitor
  auto v1{std::make_shared<xfdtd::VoltageMonitor>(
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4e-3, 0, 0},
          xfdtd::PointVector{msl_50_width, 2 * size, z}),
      xfdtd::Orientation::ZP,
      "./visualizing_data/data/quarter_wave_transformer", "v1.npy")};
  auto v2{std::make_shared<xfdtd::VoltageMonitor>(
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4.8e-3, y - 2 * size, 0},
          xfdtd::PointVector{msl_100_width, 2 * size, z}),
      xfdtd::Orientation::ZP,
      "./visualizing_data/data/quarter_wave_transformer", "v2.npy")};
  auto c1{std::make_shared<xfdtd::CurrentMonitor>(
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4e-3, 0, 0.4e-3},
          xfdtd::PointVector{9 * size, 2 * size, 1 * size}),
      xfdtd::Orientation::ZP,
      "./visualizing_data/data/quarter_wave_transformer", "c1.npy")};
  auto c2{std::make_shared<xfdtd::CurrentMonitor>(
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4.8e-3, y - 2 * size, 0.4e-3},
          xfdtd::PointVector{2 * size, 2 * size, 1 * size}),
      xfdtd::Orientation::ZP,
      "./visualizing_data/data/quarter_wave_transformer", "c2.npy")};
  auto monitors{std::vector<std::shared_ptr<Monitor>>{}};
  monitors.emplace_back(v1);
  monitors.emplace_back(v2);
  monitors.emplace_back(c1);
  monitors.emplace_back(c2);

  // define network
  auto port1{xfdtd::Port{1, 50, true, v1, c1}};
  auto port2{xfdtd::Port{2, 100, false, v2, c2}};
  auto ports{std::vector<std::unique_ptr<xfdtd::Port>>()};
  ports.emplace_back(std::make_unique<xfdtd::Port>(std::move(port1)));
  ports.emplace_back(std::make_unique<xfdtd::Port>(std::move(port2)));
  auto network{std::make_unique<xfdtd::Network>(
      std::move(ports), xt::arange<double>(2e7, 8e9, 2e7),
      "./visualizing_data/data/quarter_wave_transformer")};

  // FDTD simulation
  auto simulation{xfdtd::Simulation{
      size, std::move(objects), std::move(lumped_elements),
      std::move(boundaries), std::move(monitors), std::move(network)}};
  simulation.run(5000);
}

int main() {
  auto duration{xfdtd_example::timeSomething(quarterWaveTransformer)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration).count()};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()};
  std::cout << "It costs " << duration_in_seconds << " seconds or "
            << duration_in_milliseconds
            << " milliseconds to run Quarter Wave Transformer.\n";
}
