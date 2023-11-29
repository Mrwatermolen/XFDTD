#include "lumped_element/capacitor.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "helper.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/voltage_monitor.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/type_define.h"
#include "waveform/custom_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

void capacitor() {
  constexpr double size{1e-3};

  // create objects
  auto domain{
      xfdtd::Object{"domain",
                    std::make_unique<xfdtd::Cube>(
                        xfdtd::PointVector{-1 * size, -1 * size, -1 * size},
                        xfdtd::PointVector{6 * size, 4 * size, 4 * size}),
                    std::make_unique<xfdtd::Material>("air", 1, 1, 0, 0)}};
  auto pec_zn{xfdtd::ObjectPlane{
      "zn_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0, 0 * size},
          xfdtd::PointVector{2 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};
  auto pec_zp{xfdtd::ObjectPlane{
      "zp_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0, 1 * size},
          xfdtd::PointVector{2 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};

  // create lumped element
  constexpr double start_time{100 * 1.9065748878303884e-12};
  auto unit_step_f{[&start_time](double time) -> double {
    if (time < start_time) {
      return 0.0;
    }
    return 1.0;
  }};
  auto source{xfdtd::VoltageSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0, 0, 0},
          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::CustomWaveform>(1, unit_step_f)}};
  auto capacitor{
      xfdtd::Capacitor{std::make_unique<xfdtd::Cube>(
                           xfdtd::PointVector{1 * size, 0, 0},
                           xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
                       xfdtd::Axis::Z, 1e-11}};

  // create monitors
  auto voltage_monitor{xfdtd::VoltageMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{1 * size, 0, 0 * size},
          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
      xfdtd::Orientation::ZP, "./visualizing_data/data/capacitor/",
      "circuit_voltage.npy"}};

  auto objects{std::vector<std::shared_ptr<Object>>{}};
  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zn)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zp)));

  auto lumped_elements{std::vector<std::shared_ptr<xfdtd::LumpedElement>>()};
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::VoltageSource>(std::move(source)));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::Capacitor>(std::move(capacitor)));

  auto monitors{std::vector<std::shared_ptr<Monitor>>{}};
  monitors.emplace_back(
      std::make_shared<xfdtd::VoltageMonitor>(std::move(voltage_monitor)));

  auto simulation{xfdtd::Simulation{size, std::move(objects),
                                    std::move(lumped_elements),
                                    std::move(monitors)}};

  simulation.run(2000);
}

int main() {
  auto duration{xfdtd_example::timeSomething(capacitor)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  std::cout << "It costs " << duration_in_seconds.count() << " seconds or "
            << duration_in_milliseconds.count()
            << " milliseconds to run the Capacitor FDTD simulation."
            << "\n";
}
