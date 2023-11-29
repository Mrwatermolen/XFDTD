#include <chrono>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "helper.h"
#include "lumped_element/current_source.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/resistor.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "monitor/voltage_monitor.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/type_define.h"
#include "waveform/sinusoidal_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

void basicLumpedElementFDTD() {
  constexpr double size{1e-3};

  // create calculate domain
  auto domain{
      xfdtd::Object{"domain",
                    std::make_unique<xfdtd::Cube>(
                        xfdtd::PointVector{-3 * size, -3 * size, -3 * size},
                        xfdtd::PointVector{14 * size, 8 * size, 10 * size}),
                    std::make_unique<xfdtd::Material>("air", 1, 1, 0, 0)}};

  auto pec_plane_zn{xfdtd::ObjectPlane{
      "zn_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0 * size, 0 * size},
          xfdtd::PointVector{8 * size, 2 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};

  auto pec_plane_zp{xfdtd::ObjectPlane{
      "zp_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0 * size, 4 * size},
          xfdtd::PointVector{8 * size, 2 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};

  // create lumped elements
  auto waveform{xfdtd::SinusoidalWaveform{1, 5e8, 0}};
  auto voltage_source{xfdtd::VoltageSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0 * size, 0 * size},
          xfdtd::PointVector{1 * size, 2 * size, 4 * size}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::SinusoidalWaveform>(std::move(waveform))}};
  auto current_source{xfdtd::CurrentSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0 * size, 0 * size},
          xfdtd::PointVector{1 * size, 2 * size, 4 * size}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::SinusoidalWaveform>(std::move(waveform))}};

  auto resistance{
      xfdtd::Resistor{std::make_unique<xfdtd::Cube>(
                          xfdtd::PointVector{7 * size, 0 * size, 0 * size},
                          xfdtd::PointVector{1 * size, 2 * size, 4 * size}),
                      xfdtd::Axis::Z, 150}};

  // create monitors
  auto circuit_voltage_monitor{xfdtd::VoltageMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4 * size, 0 * size, 0 * size},
          xfdtd::PointVector{1 * size, 2 * size, 4 * size}),
      xfdtd::Orientation::ZP, "./visualizing_data/data/basic_lumped_element/",
      "circuit_voltage.npy"}};
  auto circuit_current_monitor{xfdtd::CurrentMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{4 * size, 0 * size, 3 * size},
          xfdtd::PointVector{1 * size, 2 * size, 1 * size}),
      xfdtd::Orientation::XP, "./visualizing_data/data/basic_lumped_element/",
      "circuit_current.npy"}};

  // add objects to simulation
  std::vector<std::shared_ptr<Object>> objects;
  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(domain)));
  objects.emplace_back(
      std::make_shared<xfdtd::ObjectPlane>(std::move(pec_plane_zn)));
  objects.emplace_back(
      std::make_shared<xfdtd::ObjectPlane>(std::move(pec_plane_zp)));

  // add lumped elements to simulation
  std::vector<std::shared_ptr<xfdtd::LumpedElement>> lumped_elements;
  //   lumped_elements.emplace_back(
  //       std::make_shared<xfdtd::VoltageSource>(std::move(voltage_source)));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::CurrentSource>(std::move(current_source)));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::Resistor>(std::move(resistance)));

  // add monitors to simulation
  std::vector<std::shared_ptr<Monitor>> monitors;
  monitors.emplace_back(std::make_shared<xfdtd::VoltageMonitor>(
      std::move(circuit_voltage_monitor)));
  monitors.emplace_back(std::make_shared<xfdtd::CurrentMonitor>(
      std::move(circuit_current_monitor)));

  auto simulation{xfdtd::Simulation{size, std::move(objects),
                                    std::move(lumped_elements),
                                    std::move(monitors)}};

  simulation.run(3000);
}

int main() {
  auto duration{xfdtd_example::timeSomething(basicLumpedElementFDTD)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  std::cout << "It costs " << duration_in_seconds.count() << " seconds or "
            << duration_in_milliseconds.count()
            << " milliseconds to run the Basic Lumped Element FDTD simulation."
            << "\n";
}
