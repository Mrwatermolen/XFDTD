#include "monitor/movie_monitor.h"

#include <string>

namespace xfdtd {
MovieMonitor::MovieMonitor(std::unique_ptr<Monitor> monitor,
                           size_t total_time_steps, size_t plot_step)
    : _monitor{std::move(monitor)},
      _total_time_steps{total_time_steps},
      _plot_step{plot_step} {}

MovieMonitor::MovieMonitor(const MovieMonitor &other) : Monitor(other) {
  if (&other == this) {
    return;
  }

  _monitor = other.clone();
  _total_time_steps = other._total_time_steps;
  _plot_step = other._plot_step;
}

MovieMonitor &MovieMonitor::operator=(const MovieMonitor &other) {
  if (&other == this) {
    return *this;
  }

  Monitor::operator=(other);
  _monitor = other.clone();
  _total_time_steps = other._total_time_steps;
  _plot_step = other._plot_step;
  return *this;
}

std::unique_ptr<Monitor> MovieMonitor::clone() const {
  return std::make_unique<MovieMonitor>(*this);
}

void MovieMonitor::init(
    const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
    const std::shared_ptr<const GridSpace> &grid_space,
    const std::shared_ptr<const EMF> &emf) {
  defaultInit(fdtd_basic_coff, grid_space, emf);
  _monitor->init(fdtd_basic_coff, grid_space, emf);
}

void MovieMonitor::update() {
  _monitor->update();
  auto current_time_step{getFDTDBasicCoffInstance()->getCurrentTimeStep()};
  if (current_time_step % _plot_step != 0 &&
      current_time_step != _total_time_steps - 1) {
    return;
  }
  std::stringstream ss;
  ss << std::setw(4) << std::setfill('0') << current_time_step;
  std::string name;
  ss >> name;
  _monitor->setOutputFileName(name + ".npy");
  _monitor->outputData();
}

const std::string &MovieMonitor::getOutputPath() const {
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

void MovieMonitor::outputData() {}
}  // namespace xfdtd
