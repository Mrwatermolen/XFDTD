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

#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "electromagnetic_field/electromagnetic_field.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {
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

  // std::ofstream coff_ofs{"simulation_coff_check", std::ios::out};
  // coff_ofs << "cexe max " << *std::max_element(_cexe.begin(), _cexe.end())
  //          << std::endl;
  // coff_ofs << "cexe min " << *std::min_element(_cexe.begin(), _cexe.end())
  //          << std::endl;
  // coff_ofs << "cexhz max " << *std::max_element(_cexhz.begin(), _cexhz.end())
  //          << std::endl;
  // coff_ofs << "cexhz min " << *std::min_element(_cexhz.begin(), _cexhz.end())
  //          << std::endl;
  // coff_ofs << "cexhy max " << *std::max_element(_cexhy.begin(), _cexhy.end())
  //          << std::endl;
  // coff_ofs << "cexhy min " << *std::min_element(_cexhy.begin(), _cexhy.end())
  //          << std::endl;
  // coff_ofs << "ceye max " << *std::max_element(_ceye.begin(), _ceye.end())
  //          << std::endl;
  // coff_ofs << "ceye min " << *std::min_element(_ceye.begin(), _ceye.end())
  //          << std::endl;
  // coff_ofs << "ceyhx max " << *std::max_element(_ceyhx.begin(), _ceyhx.end())
  //          << std::endl;
  // coff_ofs << "ceyhx min " << *std::min_element(_ceyhx.begin(), _ceyhx.end())
  //          << std::endl;
  // coff_ofs << "ceyhz max " << *std::max_element(_ceyhz.begin(), _ceyhz.end())
  //          << std::endl;
  // coff_ofs << "ceyhz min " << *std::min_element(_ceyhz.begin(), _ceyhz.end())
  //          << std::endl;
  // coff_ofs << "ceze max " << *std::max_element(_ceze.begin(), _ceze.end())
  //          << std::endl;
  // coff_ofs << "ceze min " << *std::min_element(_ceze.begin(), _ceze.end())
  //          << std::endl;
  // coff_ofs << "cezhx max " << *std::max_element(_cezhx.begin(), _cezhx.end())
  //          << std::endl;
  // coff_ofs << "cezhx min " << *std::min_element(_cezhx.begin(), _cezhx.end())
  //          << std::endl;
  // coff_ofs << "cezhy max " << *std::max_element(_cezhy.begin(), _cezhy.end())
  //          << std::endl;
  // coff_ofs << "cezhy min " << *std::min_element(_cezhy.begin(), _cezhy.end())
  //          << std::endl;
  // coff_ofs << "chxh max " << *std::max_element(_chxh.begin(), _chxh.end())
  //          << std::endl;
  // coff_ofs << "chxh min " << *std::min_element(_chxh.begin(), _chxh.end())
  //          << std::endl;
  // coff_ofs << "chxez max " << *std::max_element(_chxez.begin(), _chxez.end())
  //          << std::endl;
  // coff_ofs << "chxez min " << *std::min_element(_chxez.begin(), _chxez.end())
  //          << std::endl;
  // coff_ofs << "chxey max " << *std::max_element(_chxey.begin(), _chxey.end())
  //          << std::endl;
  // coff_ofs << "chxey min " << *std::min_element(_chxey.begin(), _chxey.end())
  //          << std::endl;
  // coff_ofs << "chyh max " << *std::max_element(_chyh.begin(), _chyh.end())
  //          << std::endl;
  // coff_ofs << "chyh min " << *std::min_element(_chyh.begin(), _chyh.end())
  //          << std::endl;
  // coff_ofs << "chyex max " << *std::max_element(_chyex.begin(), _chyex.end())
  //          << std::endl;
  // coff_ofs << "chyex min " << *std::min_element(_chyex.begin(), _chyex.end())
  //          << std::endl;
  // coff_ofs << "chyez max " << *std::max_element(_chyez.begin(), _chyez.end())
  //          << std::endl;
  // coff_ofs << "chyez min " << *std::min_element(_chyez.begin(), _chyez.end())
  //          << std::endl;
  // coff_ofs << "chzh max " << *std::max_element(_chzh.begin(), _chzh.end())
  //          << std::endl;
  // coff_ofs << "chzh min " << *std::min_element(_chzh.begin(), _chzh.end())
  //          << std::endl;
  // coff_ofs << "chzex max " << *std::max_element(_chzex.begin(), _chzex.end())
  //          << std::endl;
  // coff_ofs << "chzex min " << *std::min_element(_chzex.begin(), _chzex.end())
  //          << std::endl;
  // coff_ofs << "chzey max " << *std::max_element(_chzey.begin(), _chzey.end())
  //          << std::endl;
  // coff_ofs << "chzey min " << *std::min_element(_chzey.begin(), _chzey.end())
  //          << std::endl;
  // coff_ofs.close();
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
    // updateAddtiveSource();
    updateE();
    updateTFSFE();
    updateBoundaryE();
    updateTFSFIncidentField();
    updateH();
    updateTFSFH();
    updateBoundaryH();
    updateNFFFT();
    updateMonitor();
  }
  if (_nffft != nullptr) {
    _nffft->outputData();
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
  auto i{static_cast<SpatialIndex>(
      std::round((p(0) - _simulation_box->getXmin()) / _dx))};
  auto j{static_cast<SpatialIndex>(
      std::round((p(1) - _simulation_box->getYmin()) / _dy))};
  auto k{static_cast<SpatialIndex>(
      std::round((p(2) - _simulation_box->getZmin()) / _dz))};
  if (_nx == 1 && _ny == 1) {
    getEx(0, 0, k) =
        getEx(0, 0, k) +
        _cezje(0, 0, k) *
            point->getValue(_current_time_step);  // the 1d point for TFSF
  } else {
    getEz(i, j, k) =
        getEz(i, j, k) + _cezje(i, j, k) * point->getValue(_current_time_step);
  }
}

void Simulation::updateTFSFIncidentField() {
  if (_tfsf == nullptr) {
    return;
  }

  _tfsf->updateIncidentField(_current_time_step);
}

void Simulation::updateH() {
  auto& ex{getEx()};
  auto& ey{getEy()};
  auto& ez{getEz()};
  auto& hx{getHx()};
  auto& hy{getHy()};
  auto& hz{getHz()};

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        hx(i, j, k) = _chxh(i, j, k) * hx(i, j, k) +
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

  //   hx = _chxh * hx + _chxey * xt::diff(ey, 1, 2) +
  //        _chxez * xt::diff(ez, 1, 1);
  //   hy = _chyh * hy + _chyez * xt::diff(ez, 1, 0) +
  //        _chyex * xt::diff(ex, 1, 2);
  //   hz = _chzh * hz + _chzex * xt::diff(ex, 1, 1) +
  //        _chzey * xt::diff(ey, 1, 0);
}

void Simulation::updateTFSFH() {
  if (_tfsf == nullptr) {
    return;
  }

  _tfsf->updateH();
}

void Simulation::updateBoundaryH() {
  for (auto&& e : _boundaries) {
    e->updateH();
  }
}

void Simulation::updateBoundaryE() {
  for (auto&& e : _boundaries) {
    e->updateE();
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
    for (SpatialIndex k{0}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        for (SpatialIndex j{1}; j < _ny; ++j) {
          ez(i, j, k) = _ceze(i, j, k) * ez(i, j, k) +
                        _cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        _cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        }
      }
    }
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
  for (SpatialIndex j{0}; j < _ny; ++j) {
    for (SpatialIndex k{1}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        ey(i, j, k) = _ceye(i, j, k) * ey(i, j, k) +
                      _ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) +
                      _ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex k{0}; k < _nz; ++k) {
    for (SpatialIndex i{1}; i < _nx; ++i) {
      for (SpatialIndex j{1}; j < _ny; ++j) {
        ez(i, j, k) = _ceze(i, j, k) * ez(i, j, k) +
                      _cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
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

void Simulation::updateNFFFT() {
  if (_nffft == nullptr) {
    return;
  }

  _nffft->update(_current_time_step);
}

void Simulation::updateMonitor() {
  for (auto&& e : _monitors) {
    e->update(_current_time_step);
  }
}

}  // namespace xfdtd
