#include "lumped_element/inductor.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "helper.h"
#include "lumped_element/current_source.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/resistor.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/type_define.h"
#include "waveform/custom_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

void inductor() {
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

int main() {
  auto duration{xfdtd_example::timeSomething(inductor)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  std::cout << "It costs " << duration_in_seconds.count() << " seconds or "
            << duration_in_milliseconds.count()
            << " milliseconds to run the Inductor FDTD simulation."
            << "\n";
}
