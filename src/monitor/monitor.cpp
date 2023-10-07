#include "monitor/monitor.h"

#include <utility>

namespace xfdtd {
Monitor::Monitor(std::unique_ptr<Shape> shape,
                 std::filesystem::path output_dir_path,
                 std::string output_file_name)
    : _shape{std::move(shape)},
      _output_dir_path{std::move(output_dir_path)},
      _output_file_name{std::move(output_file_name)} {}

Monitor::Monitor(const Monitor& other) {
  if (&other == this) {
    return;
  }

  _shape = other._shape->clone();
  _output_dir_path = other._output_dir_path;
  _output_file_name = other._output_file_name;
  _fdtd_basic_coff = other._fdtd_basic_coff;
  _grid_space = other._grid_space;
  _emf = other._emf;
}

Monitor& Monitor::operator=(const Monitor& other) {
  if (&other == this) {
    return *this;
  }

  _shape = other._shape->clone();
  _output_dir_path = other._output_dir_path;
  _output_file_name = other._output_file_name;
  _fdtd_basic_coff = other._fdtd_basic_coff;
  _grid_space = other._grid_space;
  _emf = other._emf;
  return *this;
}

void Monitor::defaultInit(
    const std::shared_ptr<const FDTDBasicCoff>& fdtd_basic_coff,
    const std::shared_ptr<const GridSpace>& grid_space,
    const std::shared_ptr<const EMF>& emf) {
  _grid_space = grid_space;
  _fdtd_basic_coff = fdtd_basic_coff;
  _emf = emf;
}

const std::filesystem::path& Monitor::getOutputPath() const {
  return _output_dir_path;
}

const std::string& Monitor::getOutputFileName() const {
  return _output_file_name;
}

void Monitor::setOutputDirPath(const std::string& output_dir_path) {
  _output_dir_path = output_dir_path;
}

void Monitor::setOutputFileName(const std::string& output_file_name) {
  _output_file_name = output_file_name;
}

const Shape* Monitor::getShape() const { return _shape.get(); }

const GridSpace* Monitor::getGridSpaceInstance() const {
  return _grid_space.get();
}

const FDTDBasicCoff* Monitor::getFDTDBasicCoffInstance() const {
  return _fdtd_basic_coff.get();
}

const EMF* Monitor::getEMFInstance() const { return _emf.get(); }
}  // namespace xfdtd
