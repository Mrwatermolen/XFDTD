#include "monitor/monitor.h"

#include <utility>

namespace xfdtd {
Monitor::Monitor(std::unique_ptr<Shape> shape,
                 std::filesystem::path output_dir_path,
                 std::string output_file_name)
    : _shape{std::move(shape)},
      _output_dir_path{std::move(output_dir_path)},
      _output_file_name{std::move(output_file_name)} {}

const std::unique_ptr<Shape>& Monitor::getShape() const { return _shape; }

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
}  // namespace xfdtd
