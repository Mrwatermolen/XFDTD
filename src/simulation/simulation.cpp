#define XTENSOR_USE_XSIMD

#include "simulation/simulation.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>

#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "electromagnetic_field/electromagnetic_field.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {

Simulation::Simulation(double cell_size, float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _emf{std::make_shared<EMF>()} {};

Simulation::Simulation(double cell_size, ObjectArray objects,
                       BoundaryArray boundaries, std::unique_ptr<TFSF> tfsf,
                       float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _boundaries{std::move(boundaries)},
      _tfsf{std::move(tfsf)},
      _emf{std::make_shared<EMF>()} {};
Simulation::Simulation(double cell_size, ObjectArray objects,
                       BoundaryArray boundaries, std::unique_ptr<TFSF> tfsf,
                       std::unique_ptr<NFFFT> nffft, float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _boundaries{std::move(boundaries)},
      _tfsf{std::move(tfsf)},
      _nffft{std::move(nffft)},
      _emf{std::make_shared<EMF>()} {};

Simulation::Simulation(double cell_size, ObjectArray objects,
                       SourceArray sources, std::unique_ptr<TFSF> tfsf,
                       std::unique_ptr<NFFFT> nffft, BoundaryArray boundaries,
                       MonitorArray monitors, float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _sources{std::move(sources)},
      _tfsf{std::move(tfsf)},
      _nffft{std::move(nffft)},
      _boundaries{std::move(boundaries)},
      _monitors{std::move(monitors)},
      _emf{std::make_shared<EMF>()} {}

Simulation::Simulation(double cell_size, ObjectArray objects,
                       SourceArray sources, BoundaryArray boundaries,
                       MonitorArray monitors, float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _sources{std::move(sources)},
      _tfsf{nullptr},
      _boundaries{std::move(boundaries)},
      _monitors{std::move(monitors)},
      _emf{std::make_shared<EMF>()} {}

void Simulation::checkRun(size_t time_steps) {
  std::cout << "Simulation Check:" << std::endl;
  _time_steps = time_steps;
  init();
  std::ofstream ofs{"simulation_check", std::ios::out};
  std::ofstream object_ofs{"simulation_object_check", std::ios::out};
  ofs << "dt:" << _dt << " total_time_steps:" << _time_steps << std::endl;
  ofs << "dx:" << _dx << " dy:" << _dy << " dz:" << _dz << std::endl;
  ofs << "nx:" << _nx << " ny:" << _ny << " nz:" << _nz << std::endl;
  ofs << "Totoal size:" << _nx * _ny * _nz << std::endl;
  ofs << "Simulation box:" << static_cast<std::string>(*_simulation_box)
      << std::endl;
  ofs.close();

  object_ofs << "Object:" << std::endl;
  for (auto&& e : _objects) {
    object_ofs << static_cast<std::string>(*e) << std::endl;
  }
  object_ofs.close();
}

void Simulation::setCellSize(double cell_size) {
  _dx = cell_size;
  _dy = cell_size;
  _dz = cell_size;
}

void Simulation::addObject(std::unique_ptr<Object> object) {
  _objects.push_back(std::move(object));
}

void Simulation::addTFSFSource(std::unique_ptr<TFSF> tfsf) {
  _tfsf = std::move(tfsf);
}

void Simulation::addNFFFT(std::unique_ptr<NFFFT> nffft) {
  _nffft = std::move(nffft);
}

void Simulation::addMonitor(std::unique_ptr<Monitor> monitor) {
  _monitors.push_back(std::move(monitor));
}

}  // namespace xfdtd
