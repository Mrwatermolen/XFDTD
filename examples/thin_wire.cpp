#include "object/thin_wire.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "boundary/perfect_match_layer.h"
#include "helper.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "monitor/voltage_monitor.h"
#include "network/network.h"
#include "network/port.h"
#include "nffft/nffft_fd.h"
#include "object/object.h"
#include "shape/cube.h"
#include "shape/cylinder.h"
#include "simulation/simulation.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

using xfdtd::Boundary;
using xfdtd::Monitor;
using xfdtd::Object;

/**
 * @brief This is a example for simulating a half wave dipole antenna with thin
 * wire model.
 *
 */
void halfWaveDipole() {
  using namespace xfdtd;

  constexpr double cell_size{2.5e-4};
  const std::string save_path{"./visualizing_data/data/half_wave_dipole"};

  auto domain{std::make_unique<Object>(
      "domain",
      std::make_unique<Cube>(
          PointVector{-10 * cell_size, -10 * cell_size,
                      -10e-3 - 10 * cell_size},
          PointVector{20 * cell_size, 20 * cell_size, 20e-3 + 20 * cell_size}),
      Material::createAir())};

  std::vector<std::shared_ptr<Boundary>> boundaries;
  boundaries.emplace_back(std::make_unique<PML>(Orientation::XN, 8));
  boundaries.emplace_back(std::make_unique<PML>(Orientation::XP, 8));
  boundaries.emplace_back(std::make_unique<PML>(Orientation::YN, 8));
  boundaries.emplace_back(std::make_unique<PML>(Orientation::YP, 8));
  boundaries.emplace_back(std::make_unique<PML>(Orientation::ZN, 8));
  boundaries.emplace_back(std::make_unique<PML>(Orientation::ZP, 8));

  // Create half dipole antenna
  constexpr double arm_length{39 * cell_size};
  auto thin_wire_0{std::make_unique<ThinWire>(
      "arm_0", std::make_unique<Cylinder>(Axis::Z, PointVector{0, 0, 5.125e-3},
                                          5e-5, arm_length))};
  auto thin_wire_1{std::make_unique<ThinWire>(
      "arm_1", std::make_unique<Cylinder>(Axis::Z, PointVector{0, 0, -5.125e-3},
                                          5e-5, arm_length))};
  std::vector<std::shared_ptr<Object>> objects;
  objects.emplace_back(std::move(domain));
  objects.emplace_back(std::move(thin_wire_0));
  objects.emplace_back(std::move(thin_wire_1));

  // add a voltage source
  constexpr double tau{20 * cell_size / 6e8};
  constexpr double t_0{4.5 * tau};
  auto source{std::make_unique<VoltageSource>(
      std::make_unique<Cube>(PointVector{0, 0, -cell_size},
                             PointVector{cell_size, cell_size, 2 * cell_size}),
      Orientation::ZP, 50, std::make_unique<GaussianWaveform>(1, tau, t_0))};
  std::vector<std::shared_ptr<LumpedElement>> lumped_elements;
  lumped_elements.emplace_back(std::move(source));

  // add monitor to calculate S parameters
  auto v1{std::make_shared<VoltageMonitor>(
      std::make_unique<Cube>(PointVector{0, 0, -cell_size},
                             PointVector{cell_size, cell_size, 2 * cell_size}),
      Orientation::ZP, save_path, "v1.npy")};
  auto c1{std::make_shared<CurrentMonitor>(
      std::make_unique<Cube>(PointVector{0, 0, 0},
                             PointVector{cell_size, cell_size, cell_size}),
      Orientation::ZP, save_path, "c1.npy")};
  std::vector<std::shared_ptr<Monitor>> monitors;
  monitors.emplace_back(v1);
  monitors.emplace_back(c1);
  auto port_1{std::make_unique<Port>(1, 50, true, v1, c1)};
  auto ports{std::vector<std::unique_ptr<Port>>{}};
  ports.emplace_back(std::move(port_1));
  auto network{std::make_unique<Network>(
      std::move(ports), xt::arange<double>(2e7, 2e10, 2e7), save_path)};

  auto simulation{Simulation{cell_size, std::move(objects),
                             std::move(lumped_elements), std::move(boundaries),
                             std::move(monitors), std::move(network)}};

  // add a NFFFT to plot XZ plane radiation pattern
  simulation.addNFFFT(
      std::make_unique<NffftFd>(13, 13, 13, xt::xarray<double>{7e9},
                                xt::linspace<double>(0, constant::PI, 181),
                                xt::xarray<double>{0}, save_path));

  simulation.run(4000);
}

int main() {
  auto duration{xfdtd_example::timeSomething(halfWaveDipole)};
  auto duration_in_seconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  auto duration_in_milliseconds{
      std::chrono::duration_cast<std::chrono::seconds>(duration)};
  std::cout << "It costs " << duration_in_seconds.count() << " seconds or "
            << duration_in_milliseconds.count()
            << " milliseconds to run the Half Wave Dipole simulation."
            << "\n";
  return 0;
}