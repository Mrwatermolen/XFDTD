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

#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

#define DEBUG

namespace xfdtd {
Simulation::Simulation(double cell_size, ObjectArray objects,
                       SourceArray sources, std::unique_ptr<TFSF> tfsf,
                       BoundaryArray boundaries, MonitorArray monitors,
                       float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _sources{std::move(sources)},
      _tfsf{std::move(tfsf)},
      _boundaries{std::move(boundaries)},
      _monitors{std::move(monitors)} {}

void Simulation::checkRun(size_t time_steps) {
  _time_steps = time_steps;
  init();
}

void Simulation::run(size_t time_steps) {
  _time_steps = time_steps;
  init();
  std::cout << "Start running simulation..." << std::endl;
  for (size_t i{0}; i < _time_steps; ++i) {
    std::cout << "\r"
              << "Progress: " << std::setw(4) << std::setfill('0') << i + 1
              << "/" << _time_steps << std::flush;
    _current_time_step = i;
    updateAddtiveSource();
    updateTFSFIncidentField();
    updateH();
    updateTFSFH();
    updateBoundaryH();
    updateE();
    updateTFSFE();
    updateBoundaryE();
    updateMonitor();
  }
  std::cout << "\n"
            << "Simulation finished." << std::endl;
}

void Simulation::updateAddtiveSource() {
  for (auto& source : _sources) {
    handleHardPointSource(source.get());
  }
}

void Simulation::handleHardPointSource(Source* source) {
  auto point{dynamic_cast<HardPonitSource*>(source)};
  if (point == nullptr) {
    return;
  }
  auto p{point->getPoint()};
  // TODO(franzero):错误的计算
  auto k{static_cast<SpatialIndex>(
      std::round((p.z() - _simulation_box->getZmin()) / _dz))};

  getEx()(0, 0, k) =
      getEx()(0, 0, k) + _cexje(0, 0, k) * point->getValue(_current_time_step);
}

void Simulation::updateTFSFIncidentField() {
  if (_tfsf == nullptr) {
    return;
  }

  _tfsf->updateIncidentField(_current_time_step);
}

void Simulation::updateH() {
  // if (_nx == 1 && _ny == 1) {
  //   // 1D
  //   for (SpatialIndex k{0}; k < _nz; ++k) {
  //     _hy(0, 0, k) = _chyh(0, 0, k) * _hy(0, 0, k) +
  //                    _chyex(0, 0, k) * (_ex(0, 0, k + 1) - _ex(0, 0, k));
  //   }
  //   return;
  // }
  auto& ex{getEx()};
  auto& ey{getEy()};
  auto& ez{getEz()};
  auto& hx{getHx()};
  auto& hy{getHy()};
  auto& hz{getHz()};

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        getHx()(i, j, k) = _chxh(i, j, k) * getHx()(i, j, k) +
                           _chxey(i, j, k) * (ey(i, j, k + 1) - ey(i, j, k)) +
                           _chxez(i, j, k) * (ez(i, j + 1, k) - ez(i, j, k));

        hy(i, j, k) = _chyh(i, j, k) * hy(i, j, k) +
                      _chyez(i, j, k) * (ez(i + 1, j, k) - ez(i, j, k)) +
                      _chyex(i, j, k) * (ex(i, j, k + 1) - ex(i, j, k));

        hz(i, j, k) = _chzh(i, j, k) * hz(i, j, k) +
                      _chzex(i, j, k) * (ex(i, j + 1, k) - ex(i, j, k)) +
                      _chzey(i, j, k) * (ey(i + 1, j, k) - ey(i, j, k));
      }
    }
  }
}

void Simulation::updateTFSFH() {
  if (_tfsf == nullptr) {
    return;
  }

  _tfsf->updateH();
}

void Simulation::updateBoundaryH() {
  for (auto&& e : _boundaries) {
    handlePMLBoundaryH(e);
  }
}

void Simulation::updateBoundaryE() {
  for (auto&& e : _boundaries) {
    handlePMLBoundaryE(e);
  }
}

void Simulation::updateE() {
  auto& ex{getEx()};
  auto& ey{getEy()};
  auto& ez{getEz()};
  auto& hx{getHx()};
  auto& hy{getHy()};
  auto& hz{getHz()};
  // Ex的在-y,+y,-z+z方向的边界被截断
  if (_nx == 1 && _ny == 1) {
    // 1D
    for (SpatialIndex k{1}; k < _nz; ++k) {
      ex(0, 0, k) = _cexe(0, 0, k) * ex(0, 0, k) +
                    _cexhy(0, 0, k) * (hy(0, 0, k) - hy(0, 0, k - 1));
    }
    return;
  }

  if (_nz == 1) {
    // 2D
    return;
  }

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        ex(i, j, k) = _cexe(i, j, k) * ex(i, j, k) +
                      _cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      _cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      }
    }
  }
  for (SpatialIndex i{1}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        ey(i, j, k) =
            _ceye(i, j, k) * ey(i, j, k) +
            _ceyhx(i, j, k) * (getHx()(i, j, k) - getHx()(i, j, k - 1)) +
            _ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex i{1}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        ez(i, j, k) =
            _ceze(i, j, k) * ez(i, j, k) +
            _cezhx(i, j, k) * (getHx()(i, j - 1, k) - getHx()(i, j, k)) +
            _cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      }
    }
  }
}

void Simulation::updateTFSFE() {
  if (_tfsf == nullptr) {
    return;
  }

  _tfsf->updateE();
}

void Simulation::updateMonitor() {
  for (auto&& e : _monitors) {
    e->update(_current_time_step);
  }
}

void Simulation::handlePMLBoundaryH(std::shared_ptr<Boundary>& boundary) {
  auto cpml{std::dynamic_pointer_cast<PML>(boundary)};
  if (cpml == nullptr) {
    return;
  }

  auto ori{cpml->getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    cpml->updateH(getEy(), getEz(), getHy(), getHz());
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    cpml->updateH(getEz(), getEx(), getHz(), getHx());
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    cpml->updateH(getEx(), getEy(), getHx(), getHy());
  }
}

void Simulation::handlePMLBoundaryE(std::shared_ptr<Boundary>& boundary) {
  auto cpml{std::dynamic_pointer_cast<PML>(boundary)};
  if (cpml == nullptr) {
    return;
  }

  auto ori{cpml->getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    cpml->updateE(getEy(), getEz(), getHy(), getHz());
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    cpml->updateE(getEz(), getEx(), getHz(), getHx());
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    cpml->updateE(getEx(), getEy(), getHx(), getHy());
  }
}

}  // namespace xfdtd
