#include "network/network.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <xtensor/xnpy.hpp>

namespace xfdtd {

Network::Network(std::vector<std::unique_ptr<Port>> ports,
                 xt::xarray<double> frequencies, std::string output_path)
    : _ports{std::move(ports)},
      _frequencies{std::move(frequencies)},
      _output_path{std::move(output_path)} {}

void Network::init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
                   const std::shared_ptr<const GridSpace> &grid_space,
                   const std::shared_ptr<const EMF> &emf) {
  for (auto &&port : _ports) {
    port->init(fdtd_basic_coff, grid_space, emf);
  }
  for (size_t i{0}; i < _ports.size(); ++i) {
    _port_map.emplace(i, _ports[i]->getPortIndex());
  }
}

void Network::outputData() {
  if (!std::filesystem::exists(_output_path) &&
      !std::filesystem::is_directory(_output_path)) {
    try {
      std::filesystem::create_directories(_output_path);
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << _output_path << '\n';
      return;
    }
  }

  for (auto &&p : _ports) {
    p->calculateSParameter(_frequencies);
  }

  for (size_t i{0}; i < _ports.size(); ++i) {
    if (!_ports[i]->isActivatedSource()) {
      continue;
    }

    for (size_t j{0}; j < _ports.size(); ++j) {
      _s_parameters[_ports[j]->getPortIndex() * 10 +
                    _ports[i]->getPortIndex()] =
          _ports[j]->getB() / _ports[i]->getA();
    }
  }

  auto output_path{std::filesystem::path{_output_path}};
  for (const auto &e : _s_parameters) {
    auto file{output_path / ("s" + std::to_string(e.first) + ".npy")};
    auto data{xt::stack(xt::xtuple(_frequencies, e.second))};
    xt::dump_npy(file.string(), data);
  }
}

}  // namespace xfdtd
