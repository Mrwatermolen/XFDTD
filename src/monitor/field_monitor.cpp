#include "monitor/field_monitor.h"

#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <xtensor/xnpy.hpp>

#include "grid/grid_box.h"
#include "util/type_define.h"

namespace xfdtd {
TimeDomainFieldMonitor::TimeDomainFieldMonitor(
    std::unique_ptr<Cube> shape, PlaneType plane_type, EMComponent component,
    std::filesystem::path output_dir_path, std::string _output_file_name)
    : Monitor{std::move(shape), std::move(output_dir_path),
              std::move(_output_file_name)},
      _plane_type{plane_type},
      _component{component} {}

std::unique_ptr<Monitor> TimeDomainFieldMonitor::clone() const {
  return std::make_unique<TimeDomainFieldMonitor>(*this);
}

void TimeDomainFieldMonitor::init(
    const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
    const std::shared_ptr<const GridSpace> &grid_space,
    const std::shared_ptr<const EMF> &emf) {
  defaultInit(fdtd_basic_coff, grid_space, emf);
  _grid_box = *getGridSpaceInstance()->getGridBox(getShape());
  if (_plane_type == PlaneType::XY) {
    if (1 < _grid_box.getGridNumZ()) {
      throw std::runtime_error("Error: invalid grid box");
    }
  }

  if (_plane_type == PlaneType::ZX) {
    if (1 < _grid_box.getGridNumY()) {
      throw std::runtime_error("Error: invalid grid box");
    }
  }

  if (_plane_type == PlaneType::YZ) {
    if (1 < _grid_box.getGridNumX()) {
      throw std::runtime_error("Error: invalid grid box");
    }
  }
}

void TimeDomainFieldMonitor::update() {}

void TimeDomainFieldMonitor::outputData() {
  if (!std::filesystem::exists(getOutputPath()) &&
      !std::filesystem::is_directory(getOutputPath())) {
    try {
      std::filesystem::create_directory(getOutputPath());
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << getOutputPath() << '\n';
      return;
    }
  }
  const auto emf{getEMFInstance()};
  const auto &data{emf->getEMComponent(_component)};
  auto x{_grid_box.getGridOriginIndexX()};
  auto y{_grid_box.getGridOriginIndexY()};
  auto z{_grid_box.getGridOriginIndexZ()};
  auto nx{_grid_box.getGridNumX()};
  auto ny{_grid_box.getGridNumY()};
  auto nz{_grid_box.getGridNumZ()};
  if (_plane_type == PlaneType::XY) {
    xt::xarray<double> data_view{xt::zeros<double>({nx, ny})};
    for (size_t i = x; i < nx; ++i) {
      for (size_t j = y; j < ny; ++j) {
        data_view(i, j) = data(i, j, z);
      }
    }
    xt::dump_npy(getOutputPath() / getOutputFileName(), data_view);
  }
  if (_plane_type == PlaneType::ZX) {
    xt::xarray<double> data_view{xt::zeros<double>({nx, nz})};
    for (size_t k = z; k < nz; ++k) {
      for (size_t i = x; i < nx; ++i) {
        data_view(i, k) = data(i, y, k);
      }
    }
    xt::dump_npy(getOutputPath() / getOutputFileName(), data_view);
  }
  if (_plane_type == PlaneType::YZ) {
    xt::xarray<double> data_view{xt::zeros<double>({ny, nz})};
    for (size_t j = y; j < ny; ++j) {
      for (size_t k = z; k < nz; ++k) {
        data_view(j, k) = data(x, j, k);
      }
    }
    xt::dump_npy(getOutputPath() / getOutputFileName(), data_view);
  }
}

}  // namespace xfdtd
