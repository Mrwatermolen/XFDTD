#include "simulation/simulation.h"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <xtensor/xadapt.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xview.hpp>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"

namespace xfdtd {

Simulation::Simulation(double cell_size, float cfl)
    : _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {};

Simulation::Simulation(double cell_size,
                       std::vector<std::shared_ptr<Object>> objects,
                       std::vector<std::shared_ptr<Boundary>> boundaries,
                       std::unique_ptr<TFSF> tfsf, float cfl)
    : _objects{std::move(objects)},
      _boundaries{std::move(boundaries)},
      _tfsf{std::move(tfsf)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {};

Simulation::Simulation(double cell_size,
                       std::vector<std::shared_ptr<Object>> objects,
                       std::vector<std::shared_ptr<Boundary>> boundaries,
                       std::unique_ptr<TFSF> tfsf, std::unique_ptr<NFFFT> nffft,
                       float cfl)
    : _objects{std::move(objects)},
      _boundaries{std::move(boundaries)},
      _tfsf{std::move(tfsf)},
      _nffft{std::move(nffft)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {};

Simulation::Simulation(
    double cell_size, std::vector<std::shared_ptr<Object>> objects,
    std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
    std::vector<std::shared_ptr<Boundary>> boundaries,
    std::vector<std::shared_ptr<Monitor>> monitors, float cfl)
    : _objects{std::move(objects)},
      _lumped_elements{std::move(lumped_elements)},
      _boundaries{std::move(boundaries)},
      _monitors{std::move(monitors)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {};

Simulation::Simulation(double cell_size,
                       std::vector<std::shared_ptr<Object>> objects,
                       std::vector<std::shared_ptr<Monitor>> monitors,
                       float cfl)
    : _objects{std::move(objects)},
      _monitors{std::move(monitors)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {};

Simulation::Simulation(
    double cell_size, std::vector<std::shared_ptr<Object>> objects,
    std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
    std::vector<std::shared_ptr<Monitor>> monitors, float cfl)
    : _objects{std::move(objects)},
      _lumped_elements{std::move(lumped_elements)},
      _monitors{std::move(monitors)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {}

Simulation::Simulation(
    double cell_size, std::vector<std::shared_ptr<Object>> objects,
    std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
    std::vector<std::shared_ptr<Boundary>> boundaries,
    std::vector<std::shared_ptr<Monitor>> monitors,
    std::unique_ptr<Network> network, float cfl)
    : _objects{std::move(objects)},
      _lumped_elements{std::move(lumped_elements)},
      _boundaries{std::move(boundaries)},
      _monitors{std::move(monitors)},
      _network{std::move(network)},
      _emf{std::make_shared<EMF>()},
      _grid_space{std::make_shared<GridSpace>(cell_size, cell_size, cell_size)},
      _fdtd_basic_coff{std::make_shared<FDTDBasicCoff>(cfl)} {}

void Simulation::checkRun(size_t total_time_steps) {
  std::cout << "Simulation Check:" << '\n';
  _total_time_steps = total_time_steps;
  init();
  std::ofstream ofs{"simulation_check", std::ios::out};
  std::ofstream object_ofs{"simulation_object_check", std::ios::out};
  ofs << "nx:" << _nx << " ny:" << _ny << " nz:" << _nz << '\n';
  ofs << "Total size:" << _nx * _ny * _nz << '\n';
  ofs.close();

  object_ofs << "Object:" << '\n';
  for (auto&& e : _objects) {
    object_ofs << static_cast<std::string>(*e) << '\n';
  }
  object_ofs.close();
}

void Simulation::addObject(std::shared_ptr<Object> object) {
  _objects.push_back(std::move(object));
}

void Simulation::addBoundary(std::shared_ptr<Boundary> boundary) {
  _boundaries.push_back(std::move(boundary));
}

void Simulation::addLumpedElement(
    std::shared_ptr<LumpedElement> lumped_element) {
  _lumped_elements.push_back(std::move(lumped_element));
}

void Simulation::addTFSFSource(std::shared_ptr<TFSF> tfsf) {
  _tfsf = std::move(tfsf);
}

void Simulation::addNFFFT(std::shared_ptr<NFFFT> nffft) {
  _nffft = std::move(nffft);
}

void Simulation::addMonitor(std::shared_ptr<Monitor> monitor) {
  _monitors.push_back(std::move(monitor));
}

void Simulation::addNetwork(std::shared_ptr<Network> network) {
  _network = std::move(network);
}

void Simulation::outputData() {
  std::cout << "Simulation starts to output data...\n";
  if (_nffft != nullptr) {
    _nffft->outputData();
  }
  for (auto&& e : _monitors) {
    e->outputData();
  }
  if (_network != nullptr) {
    _network->outputData();
  }
  // xt::dump_npy("./visualizing_data/data/cexe", _fdtd_basic_coff->getCexe());
  // xt::dump_npy("./visualizing_data/data/cexhy",
  // _fdtd_basic_coff->getCexhy());
  // xt::dump_npy("./visualizing_data/data/cexhz",
  // _fdtd_basic_coff->getCexhz()); xt::dump_npy("./visualizing_data/data/ceye",
  // _fdtd_basic_coff->getCeye()); xt::dump_npy("./visualizing_data/data/ceyhz",
  // _fdtd_basic_coff->getCeyhz());
  // xt::dump_npy("./visualizing_data/data/ceyhx",
  // _fdtd_basic_coff->getCeyhx()); xt::dump_npy("./visualizing_data/data/ceze",
  // _fdtd_basic_coff->getCeze()); xt::dump_npy("./visualizing_data/data/cezhx",
  // _fdtd_basic_coff->getCezhx());
  // xt::dump_npy("./visualizing_data/data/cezhy",
  // _fdtd_basic_coff->getCezhy());
}

void Simulation::outputTFSFIncidentWaveFastFourierTransform(
    const std::filesystem::path& path) {
  if (_tfsf == nullptr) {
    return;
  }
  if (!std::filesystem::exists(path) || !std::filesystem::is_directory(path)) {
    try {
      std::filesystem::create_directories(path);
    } catch (std::filesystem::filesystem_error& e) {
      std::cerr << "Error: " << e.what() << '\n';
      return;
    }
  }
  auto f_file{(path / "frequencies.npy").string()};
  auto v_file{(path / "value.npy").string()};
  auto [frequencies, v] = _tfsf->getIncidentWaveFastFourierTransform();
  xt::dump_npy(f_file, frequencies);
  xt::dump_npy(v_file, v);
}

}  // namespace xfdtd
