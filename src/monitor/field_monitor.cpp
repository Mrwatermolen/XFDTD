#include "monitor/field_monitor.h"

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include "electromagnetic.h"
#include "simulation/yee_cell.h"
#include "util/type_define.h"

namespace xfdtd {
TimeDomainFieldMonitor::TimeDomainFieldMonitor(
    std::unique_ptr<Shape> shape, PlaneType plane_type, EMComponent component,
    std::filesystem::path output_dir_path, std::string _output_file_name)
    : Monitor{std::move(shape), std::move(output_dir_path),
              std::move(_output_file_name)},
      _plane_type{plane_type},
      _component{component} {}

void TimeDomainFieldMonitor::setYeeCells(const YeeCellArray& yee_cells) {
  _yee_cells = yee_cells;
}

void TimeDomainFieldMonitor::setYeeCells(YeeCellArray&& yee_cells) {
  _yee_cells = std::move(yee_cells);
}

void TimeDomainFieldMonitor::update(size_t current_time_step) {}

void TimeDomainFieldMonitor::outputData() {
  if (!std::filesystem::exists(getOutputPath()) &&
      !std::filesystem::is_directory(getOutputPath())) {
    try {
      std::filesystem::create_directory(getOutputPath());
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << getOutputPath()
                << std::endl;
      return;
    }
  }
  std::ofstream output_file{getOutputPath() / getOutputFileName()};
  if (!output_file.is_open()) {
    std::cerr << "Error: cannot open file " << getOutputFileName() << std::endl;
    return;
  }
  output_file << "Field Type:" << &_component;

  if (_plane_type == PlaneType::XY) {
    SpatialIndex counter{-1};
    for (const auto& e : _yee_cells) {
      auto [x, y, z] = e->getGridXYZIndex();
      if (counter != x) {
        output_file << std::endl;
        counter = x;
      }
      output_file << getEMComponent(_component, x, y, z) << " ";
    }
  } else if (_plane_type == PlaneType::YZ) {
    SpatialIndex counter{-1};
    for (const auto& e : _yee_cells) {
      auto [x, y, z] = e->getGridXYZIndex();
      if (counter != y) {
        output_file << std::endl;
        counter = y;
      }
      output_file << getEMComponent(_component, x, y, z) << " ";
    }
  } else if (_plane_type == PlaneType::ZX) {
    SpatialIndex counter{-1};
    for (const auto& e : _yee_cells) {
      auto [x, y, z] = e->getGridXYZIndex();
      if (counter != x) {
        output_file << std::endl;
        counter = x;
      }
      output_file << getEMComponent(_component, x, y, z) << " ";
    }
  }
  output_file.close();
}

}  // namespace xfdtd
