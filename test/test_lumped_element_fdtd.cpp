/**
 * @file test_lumped_element_fdtd.cpp
 * @brief Example for simple lumped element FDTD simulation.Here we use pec
 * plane to simulate the wire.
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <chrono>
#include <exception>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "boundary/boundary.h"
#include "lumped_element/capacitor.h"
#include "lumped_element/current_source.h"
#include "lumped_element/inductor.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/resistor.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "monitor/monitor.h"
#include "monitor/voltage_monitor.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/type_define.h"
#include "waveform/cosine_modulated_gaussian_waveform.h"
#include "waveform/custom_waveform.h"
#include "waveform/sinusoidal_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

void testBasicLumpedElementFDTD() {
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

void testCapacitor() {
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

void testInductor() {
  constexpr double size{1e-3};

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
          xfdtd::PointVector{4 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};
  auto pec_zp{xfdtd::ObjectPlane{
      "zp_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0, 1 * size},
          xfdtd::PointVector{4 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};

  constexpr double start_time{100 * 1.9065748878303884e-12};
  auto unit_step_f{[&start_time](double time) -> double {
    if (time < start_time) {
      return 0.0;
    }
    return 1.0;
  }};
  auto current_source{xfdtd::CurrentSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0, 0, 0},
          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::CustomWaveform>(1, unit_step_f)}};
  auto inductor{
      xfdtd::Inductor{std::make_unique<xfdtd::Cube>(
                          xfdtd::PointVector{1 * size, 0, 0},
                          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
                      xfdtd::Axis::Z, 1e-8}};
  auto resistance{
      xfdtd::Resistor{std::make_unique<xfdtd::Cube>(
                          xfdtd::PointVector{1 * size, 0, 0 * size},
                          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
                      xfdtd::Axis::Z, 50}};

  auto current_monitor{xfdtd::CurrentMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0 * size, 1 * size},
          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
      xfdtd::Orientation::XP, "./visualizing_data/data/inductor/",
      "current.npy"}};

  auto objects{std::vector<std::shared_ptr<Object>>{}};
  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zn)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zp)));

  auto lumped_elements{std::vector<std::shared_ptr<xfdtd::LumpedElement>>()};
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::CurrentSource>(std::move(current_source)));
  //   lumped_elements.emplace_back(std::make_shared<xfdtd::Resistor>(resistance));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::Inductor>(std::move(inductor)));

  auto monitors{std::vector<std::shared_ptr<Monitor>>{}};
  monitors.emplace_back(
      std::make_shared<xfdtd::CurrentMonitor>(std::move(current_monitor)));

  auto simulation{xfdtd::Simulation{size, std::move(objects),
                                    std::move(lumped_elements),
                                    std::move(monitors)}};

  simulation.run(2000);
}

void testRLCCircuitFDTD() {
  constexpr double size{1e-3};

  // create calculate domain
  auto domain{
      xfdtd::Object{"domain",
                    std::make_unique<xfdtd::Cube>(
                        xfdtd::PointVector{-1 * size, -1 * size, -1 * size},
                        xfdtd::PointVector{6 * size, 4 * size, 4 * size}),
                    std::make_unique<xfdtd::Material>("air", 1, 1, 0, 0)}

  };
  auto pec_zn{xfdtd::ObjectPlane{
      "zn_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0, 0 * size},
          xfdtd::PointVector{4 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};
  auto pec_zp{xfdtd::ObjectPlane{
      "zp_plane",
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0 * size, 0, 2 * size},
          xfdtd::PointVector{4 * size, 1 * size, 0 * size}),
      std::make_unique<xfdtd::Material>("pec", 1, 1, 1e10, 0)}};

  // create lumped elements
  constexpr double bandwidth{4e9};
  constexpr double tau{0.996 / bandwidth};
  constexpr double t_0{4.5 * tau};
  auto voltage_source{xfdtd::VoltageSource{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{0, 0, 0},
          xfdtd::PointVector{1 * size, 1 * size, 2 * size}),
      xfdtd::Orientation::ZP, 50,
      std::make_unique<xfdtd::CosineModulatedGaussianWaveform>(1, tau, t_0,
                                                               2e9)}};

  auto inductor{
      xfdtd::Inductor{std::make_unique<xfdtd::Cube>(
                          xfdtd::PointVector{3 * size, 0, 0},
                          xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
                      xfdtd::Axis::Z, 1e-8}};

  auto capacitor{
      xfdtd::Capacitor{std::make_unique<xfdtd::Cube>(
                           xfdtd::PointVector{3 * size, 0, 1 * size},
                           xfdtd::PointVector{1 * size, 1 * size, 1 * size}),
                       xfdtd::Axis::Z, 1e-11}};

  // create monitors
  auto circuit_voltage_monitor{xfdtd::VoltageMonitor{
      std::make_unique<xfdtd::Cube>(
          xfdtd::PointVector{3 * size, 0, 0 * size},
          xfdtd::PointVector{1 * size, 1 * size, 2 * size}),
      xfdtd::Orientation::ZP, "./visualizing_data/data/rlc_circuit/",
      "circuit_voltage.npy"}};

  auto objects{std::vector<std::shared_ptr<Object>>{}};
  objects.emplace_back(std::make_shared<xfdtd::Object>(std::move(domain)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zn)));
  objects.emplace_back(std::make_shared<xfdtd::ObjectPlane>(std::move(pec_zp)));

  auto lumped_elements{std::vector<std::shared_ptr<xfdtd::LumpedElement>>()};
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::VoltageSource>(std::move(voltage_source)));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::Inductor>(std::move(inductor)));
  lumped_elements.emplace_back(
      std::make_shared<xfdtd::Capacitor>(std::move(capacitor)));

  auto monitors{std::vector<std::shared_ptr<Monitor>>{}};
  monitors.emplace_back(std::make_shared<xfdtd::VoltageMonitor>(
      std::move(circuit_voltage_monitor)));

  auto simulation{xfdtd::Simulation{size, std::move(objects),
                                    std::move(lumped_elements),
                                    std::move(monitors)}};
  simulation.run(2000);
}

int main() {
  auto start_time{std::chrono::steady_clock::now()};
  try {
    testRLCCircuitFDTD();
  } catch (std::exception e) {
    std::cerr << "Exception: " << e.what() << '\n';
  }
  auto end_time{std::chrono::steady_clock::now()};
  auto elapsed_time{end_time - start_time};
  auto ms{std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time)
              .count()};
  auto s{
      std::chrono::duration_cast<std::chrono::seconds>(elapsed_time).count()};
  std::cout << "Elapsed time: " << ms << " ms" << '\n';
  std::cout << "Elapsed time: " << s << " s" << '\n';
}
