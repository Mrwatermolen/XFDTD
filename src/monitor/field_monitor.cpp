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
TimeDomainFieldMonitor::TimeDomainFieldMonitor(std::unique_ptr<Cube> shape,
                                               PlaneType plane_type,
                                               EMComponent component,
                                               std::string output_dir_path,
                                               std::string output_file_name)
    : Monitor{std::move(shape), std::move(output_dir_path),
              std::move(output_file_name)},
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
    _z =
        getGridSpaceInstance()->getGridCenterZ(_grid_box.getGridOriginIndexZ());
    _x = xt::make_lambda_xfunction(
        [&grid_space](const auto &e) { return grid_space->getGridCenterX(e); },
        xt::arange<size_t>(_grid_box.getGridOriginIndexX(),
                           _grid_box.getGridEndIndexX()));
    _y = xt::make_lambda_xfunction(
        [&grid_space](const auto &e) { return grid_space->getGridCenterY(e); },
        xt::arange<size_t>(_grid_box.getGridOriginIndexY(),
                           _grid_box.getGridEndIndexY()));
  }

  if (_plane_type == PlaneType::ZX) {
    if (1 < _grid_box.getGridNumY()) {
      throw std::runtime_error("Error: invalid grid box");

      _y = getGridSpaceInstance()->getGridCenterY(
          _grid_box.getGridOriginIndexY());
      _z = xt::make_lambda_xfunction(
          [&grid_space](const auto &e) {
            return grid_space->getGridCenterZ(e);
          },
          xt::arange<size_t>(_grid_box.getGridOriginIndexZ(),
                             _grid_box.getGridEndIndexZ()));
      _x = xt::make_lambda_xfunction(
          [&grid_space](const auto &e) {
            return grid_space->getGridCenterX(e);
          },
          xt::arange<size_t>(_grid_box.getGridOriginIndexX(),
                             _grid_box.getGridEndIndexX()));
    }
  }

  if (_plane_type == PlaneType::YZ) {
    if (1 < _grid_box.getGridNumX()) {
      throw std::runtime_error("Error: invalid grid box");
    }
    _x =
        getGridSpaceInstance()->getGridCenterX(_grid_box.getGridOriginIndexX());
    _y = xt::make_lambda_xfunction(
        [&grid_space](const auto &e) { return grid_space->getGridCenterY(e); },
        xt::arange<size_t>(_grid_box.getGridOriginIndexY(),
                           _grid_box.getGridEndIndexY()));
    _z = xt::make_lambda_xfunction(
        [&grid_space](const auto &e) { return grid_space->getGridCenterZ(e); },
        xt::arange<size_t>(_grid_box.getGridOriginIndexZ(),
                           _grid_box.getGridEndIndexZ()));
  }
}

void TimeDomainFieldMonitor::update() {}

void TimeDomainFieldMonitor::outputData() {
  auto output_path{std::filesystem::path{getOutputPath()}};
  if (!std::filesystem::exists(output_path) &&
      !std::filesystem::is_directory(output_path)) {
    try {
      std::filesystem::create_directories(output_path);
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << output_path.string()
                << '\n';
      return;
    }
  }

  const auto emf{getEMFInstance()};
  const auto &data{emf->getEMComponent(_component)};
  auto is{_grid_box.getGridOriginIndexX()};
  auto js{_grid_box.getGridOriginIndexY()};
  auto ks{_grid_box.getGridOriginIndexZ()};
  auto nx{_grid_box.getGridNumX()};
  auto ny{_grid_box.getGridNumY()};
  auto nz{_grid_box.getGridNumZ()};
  nx = nx == 0 ? 1 : nx;
  ny = ny == 0 ? 1 : ny;
  nz = nz == 0 ? 1 : nz;
  auto x_range{xt::range(is, is + nx)};
  auto y_range{xt::range(js, js + ny)};
  auto z_range{xt::range(ks, ks + nz)};
  auto data_view{xt::view(data, x_range, y_range, z_range)};
  xt::dump_npy(
      std::filesystem::path(output_path / getOutputFileName()).string(),
      data_view);
}

}  // namespace xfdtd
