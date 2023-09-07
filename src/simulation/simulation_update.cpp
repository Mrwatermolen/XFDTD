#include "additive_source/hard_point_source.h"
#include "additive_source/source.h"
#include "simulation/simulation.h"
#include "simulation/yee_cell.h"

namespace xfdtd {
void Simulation::run(size_t time_steps) {
  _time_steps = time_steps;
  init();
  std::cout << "Start running simulation..." << std::endl;
  for (size_t i{0}; i < _time_steps; ++i) {
    std::cout << "\r"
              << "Progress: " << std::setw(4) << std::setfill('0') << i + 1
              << "/" << _time_steps << std::flush;
    _current_time_step = i;
    if (!_is_exist_dispersive_material) {
      updateE();
    } else {
      updateEWithDispersiveMaterial();
    }
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

void Simulation::updateEWithDispersiveMaterial() {
  auto& ex{getEx()};
  auto& ey{getEy()};
  auto& ez{getEz()};
  auto& hx{getHx()};
  auto& hy{getHy()};
  auto& hz{getHz()};

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        auto material_index{
            _grid_space[(i * (_ny * _nz) + j * _nz + k)]->getMaterialIndex()};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEx(i, j, k);
          continue;
        }
        ex(i, j, k) = _cexe(i, j, k) * ex(i, j, k) +
                      _cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      _cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      }
    }
  }
  for (SpatialIndex j{0}; j < _ny; ++j) {
    for (SpatialIndex k{1}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        auto material_index{
            _grid_space[(i * (_ny * _nz) + j * _nz + k)]->getMaterialIndex()};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEy(i, j, k);
          continue;
        }
        ey(i, j, k) = _ceye(i, j, k) * ey(i, j, k) +
                      _ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) +
                      _ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex k{0}; k < _nz; ++k) {
    for (SpatialIndex i{1}; i < _nx; ++i) {
      for (SpatialIndex j{1}; j < _ny; ++j) {
        auto material_index{
            _grid_space[(i * (_ny * _nz) + j * _nz + k)]->getMaterialIndex()};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEz(i, j, k);
          continue;
        }
        ez(i, j, k) = _ceze(i, j, k) * ez(i, j, k) +
                      _cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      _cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      }
    }
  }
}
}  // namespace xfdtd