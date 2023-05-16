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

#define DEBUG

namespace xfdtd {
Simulation::Simulation(double cell_size, ObjectArray objects,
                       SourceArray sources, BoundaryArray boundaries, float cfl)
    : _dx{cell_size},
      _dy{cell_size},
      _dz{cell_size},
      _cfl{cfl},
      _objects{std::move(objects)},
      _sources{std::move(sources)},
      _boundaries{std::move(boundaries)} {}

void Simulation::setMonitor() {}

void Simulation::checkRun(size_t time_steps) {
  _time_steps = time_steps;
  init();
}

void Simulation::run(size_t time_steps) {
  _time_steps = time_steps;
  init();
  std::cout << "Start running simulation..." << std::endl;
  auto oup_dir = std::filesystem::absolute("output");
  std::cout << "Out dir: " << oup_dir.c_str() << std::endl;
  if (!std::filesystem::exists(oup_dir) &&
      !std::filesystem::is_directory(oup_dir)) {
    try {
      std::filesystem::create_directory(oup_dir);
    } catch (std::exception e) {
      std::cerr << e.what() << std::endl;
    }
  }
  std::cout << "\n";
  for (size_t i{0}; i < _time_steps; ++i) {
    std::cout << "\r"
              << "Progress: " << std::setw(4) << std::setfill('0') << i + 1
              << "/" << _time_steps << std::flush;
    _current_time_step = i;
    updateAddtiveSource();
    updateH();
    updateBoundaryH();
    updateE();
    updateBoundaryE();
#ifdef DEBUG
    if (i % 30 == 0 || i == time_steps - 1) {
      std::stringstream ss;
      ss << std::setw(4) << std::setfill('0') << i << "\n";
      std::string s;
      ss >> s;
      std::ofstream f(std::filesystem::absolute(oup_dir) /
                      ("Ex-" + s + ".dat"));
      if (!f.is_open()) {
        return;
      }
      for (auto k{0}; k < _nz; ++k) {
        f << k << "\t" << _ex(0, 0, k) << "\n";
      }
      f.close();
    }
#endif  // DEBUG
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

  _ex(0, 0, k) =
      _ex(0, 0, k) + _cexje(0, 0, k) * point->getValue(_current_time_step);
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
  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        _hx(i, j, k) = _chxh(i, j, k) * _hx(i, j, k) +
                       _chxey(i, j, k) * (_ey(i, j, k + 1) - _ey(i, j, k)) +
                       _chxez(i, j, k) * (_ez(i, j + 1, k) - _ez(i, j, k));

        _hy(i, j, k) = _chyh(i, j, k) * _hy(i, j, k) +
                       _chyez(i, j, k) * (_ez(i + 1, j, k) - _ez(i, j, k)) +
                       _chyex(i, j, k) * (_ex(i, j, k + 1) - _ex(i, j, k));

        _hz(i, j, k) = _chzh(i, j, k) * _hz(i, j, k) +
                       _chzex(i, j, k) * (_ex(i, j + 1, k) - _ex(i, j, k)) +
                       _chzey(i, j, k) * (_ey(i + 1, j, k) - _ey(i, j, k));
      }
    }
  }
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
  // Ex的在-y,+y,-z+z方向的边界被截断
  if (_nx == 1 && _ny == 1) {
    // 1D
    for (SpatialIndex k{1}; k < _nz; ++k) {
      _ex(0, 0, k) = _cexe(0, 0, k) * _ex(0, 0, k) +
                     _cexhy(0, 0, k) * (_hy(0, 0, k) - _hy(0, 0, k - 1));
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
        _ex(i, j, k) = _cexe(i, j, k) * _ex(i, j, k) +
                       _cexhy(i, j, k) * (_hy(i, j, k) - _hy(i, j, k - 1)) +
                       _cexhz(i, j, k) * (_hz(i, j, k) - _hz(i, j - 1, k));
      }
    }
  }
  for (SpatialIndex i{1}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        _ey(i, j, k) = _ceye(i, j, k) * _ey(i, j, k) +
                       _ceyhx(i, j, k) * (_hx(i, j, k) - _hx(i, j, k - 1)) +
                       _ceyhz(i, j, k) * (_hz(i, j, k) - _hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex i{1}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        _ez(i, j, k) = _ceze(i, j, k) * _ez(i, j, k) +
                       _cezhx(i, j, k) * (_hx(i, j - 1, k) - _hx(i, j, k)) +
                       _cezhy(i, j, k) * (_hy(i, j, k) - _hy(i - 1, j, k));
      }
    }
  }
}

void Simulation::output() {
  std::ofstream data_file("1d_fdtd.dat");
  data_file.close();
}

const EFTA& Simulation::getEx() const { return _ex; }
const EFTA& Simulation::getEy() const { return _ey; }
const EFTA& Simulation::getEz() const { return _ez; }
const EFTA& Simulation::getHx() const { return _hx; }
const EFTA& Simulation::getHy() const { return _hy; }
const EFTA& Simulation::getHz() const { return _hz; }

void Simulation::handlePMLBoundaryH(std::shared_ptr<Boundary>& boundary) {
  auto cpml{std::dynamic_pointer_cast<PML>(boundary)};
  if (cpml == nullptr) {
    return;
  }

  auto ori{cpml->getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    cpml->updateH(_ey, _ez, _hy, _hz);
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    cpml->updateH(_ez, _ex, _hz, _hx);
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    cpml->updateH(_ex, _ey, _hx, _hy);
  }
}

void Simulation::handlePMLBoundaryE(std::shared_ptr<Boundary>& boundary) {
  auto cpml{std::dynamic_pointer_cast<PML>(boundary)};
  if (cpml == nullptr) {
    return;
  }

  auto ori{cpml->getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    cpml->updateE(_ey, _ez, _hy, _hz);
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    cpml->updateE(_ez, _ex, _hz, _hx);
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    cpml->updateE(_ex, _ey, _hx, _hy);
  }
}

}  // namespace xfdtd
