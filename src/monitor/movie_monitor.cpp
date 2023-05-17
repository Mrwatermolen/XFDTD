#include "monitor/movie_monitor.h"

#include <iomanip>
#include <sstream>
#include <utility>

#include "util/type_define.h"

namespace xfdtd {
MovieMonitor::MovieMonitor(std::unique_ptr<Monitor> monitor,
                           size_t total_time_steps, size_t plot_step)
    : _monitor{std::move(monitor)},
      _total_time_steps{total_time_steps},
      _plot_step{plot_step} {}

const std::unique_ptr<Shape> &MovieMonitor::getShape() const {
  return _monitor->getShape();
}

const std::filesystem::path &MovieMonitor::getOutputPath() const {
  return _monitor->getOutputPath();
}
const std::string &MovieMonitor::getOutputFileName() const {
  return _monitor->getOutputFileName();
}

void MovieMonitor::setOutputDirPath(const std::string &output_dir_path) {
  _monitor->setOutputDirPath(output_dir_path);
}

void MovieMonitor::setOutputFileName(const std::string &output_file_name) {
  _monitor->setOutputFileName(output_file_name);
}

void MovieMonitor::setYeeCells(const YeeCellArray &yee_cells) {
  _monitor->setYeeCells(yee_cells);
}
void MovieMonitor::setYeeCells(YeeCellArray &&yee_cells) {
  _monitor->setYeeCells(std::move(yee_cells));
}

void MovieMonitor::update(size_t current_time_step) {
  if (current_time_step % _plot_step == 0 ||
      current_time_step == _total_time_steps - 1) {
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << current_time_step;
    std::string name;
    ss >> name;
    _monitor->update(current_time_step);
    _monitor->setOutputFileName(name + ".dat");
    _monitor->outputData();
  }
}

void MovieMonitor::outputData() {}
}  // namespace xfdtd