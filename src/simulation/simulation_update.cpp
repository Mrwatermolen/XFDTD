#include <xtensor/xnpy.hpp>

#include "simulation/simulation.h"
namespace xfdtd {
void Simulation::run(size_t total_time_steps) {
  _total_time_steps = total_time_steps;
  init();
  std::cout << "Start running simulation..." << '\n';
  for (size_t i{_current_time_step}; i < _total_time_steps; ++i) {
    std::cout << "\r"
              << "Progress: " << std::setw(4) << std::setfill('0') << i + 1
              << "/" << _total_time_steps << std::flush;
    _current_time_step = i;
    _fdtd_basic_coff->setCurrentTimeStep(_current_time_step);
    // IMPORTANT: The order of update E/H field is important. It will affect the
    // result. You will get totally different result if you change the order.

    // Time: (current_time_step) * dt
    if (!_is_exist_dispersive_material) {
      updateE();
    } else {
      updateEWithDispersiveMaterial();
    }
    correctE();
    updateBoundaryE();

    // Time: (current_time_step + 0.5) * dt
    updateH();
    correctH();
    updateBoundaryH();
    updateTFSFIncidentField();

    // finish calculation. Record
    updateNFFFT();
    updateMonitor();
  }
  std::cout << "\n";
  outputData();
  std::cout << "Simulation finished." << '\n';
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
  const auto& chxh{getChxh()};
  const auto& chxey{getChxey()};
  const auto& chxez{getChxez()};
  const auto& chyh{getChyh()};
  const auto& chyez{getChyez()};
  const auto& chyex{getChyex()};
  const auto& chzh{getChzh()};
  const auto& chzex{getChzex()};
  const auto& chzey{getChzey()};
  // for (SpatialIndex i{0}; i < _nx; ++i) {
  //   for (SpatialIndex j{0}; j < _ny; ++j) {
  //     for (SpatialIndex k{0}; k < _nz; ++k) {
  //       hx.at(i, j, k) =
  //           chxh.at(i, j, k) * hx.at(i, j, k) +
  //           chxey.at(i, j, k) * (ey.at(i, j, k + 1) - ey.at(i, j, k)) +
  //           chxez.at(i, j, k) * (ez.at(i, j + 1, k) - ez.at(i, j, k));
  //       hy.at(i, j, k) =
  //           chyh.at(i, j, k) * hy.at(i, j, k) +
  //           chyez.at(i, j, k) * (ez.at(i + 1, j, k) - ez.at(i, j, k)) +
  //           chyex.at(i, j, k) * (ex.at(i, j, k + 1) - ex.at(i, j, k));
  //       hz.at(i, j, k) =
  //           chzh.at(i, j, k) * hz.at(i, j, k) +
  //           chzex.at(i, j, k) * (ex.at(i, j + 1, k) - ex.at(i, j, k)) +
  //           chzey.at(i, j, k) * (ey.at(i + 1, j, k) - ey.at(i, j, k));
  //     }
  //   }
  // }
  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        hx(i, j, k) = chxh(i, j, k) * hx(i, j, k) +
                      chxey(i, j, k) * (ey(i, j, k + 1) - ey(i, j, k)) +
                      chxez(i, j, k) * (ez(i, j + 1, k) - ez(i, j, k));

        hy(i, j, k) = chyh(i, j, k) * hy(i, j, k) +
                      chyez(i, j, k) * (ez(i + 1, j, k) - ez(i, j, k)) +
                      chyex(i, j, k) * (ex(i, j, k + 1) - ex(i, j, k));

        hz(i, j, k) = chzh(i, j, k) * hz(i, j, k) +
                      chzex(i, j, k) * (ex(i, j + 1, k) - ex(i, j, k)) +
                      chzey(i, j, k) * (ey(i + 1, j, k) - ey(i, j, k));
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
  const auto& cexe{getCexe()};
  const auto& cexhy{getCexhy()};
  const auto& cexhz{getCexhz()};
  const auto& ceye{getCeye()};
  const auto& ceyhz{getCeyhz()};
  const auto& ceyhx{getCeyhx()};
  const auto& ceze{getCeze()};
  const auto& cezhx{getCezhx()};
  const auto& cezhy{getCezhy()};

  // Ex的在-y,+y,-z+z方向的边界被截断
  if (_nx == 1 && _ny == 1) {
    // 1D
    for (SpatialIndex k{1}; k < _nz; ++k) {
      ex(0, 0, k) = cexe(0, 0, k) * ex(0, 0, k) +
                    cexhy(0, 0, k) * (hy(0, 0, k) - hy(0, 0, k - 1));
    }
    return;
  }

  if (_nz == 1) {
    for (SpatialIndex k{0}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        for (SpatialIndex j{1}; j < _ny; ++j) {
          ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                        cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                        cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
        }
      }
    }
    return;
  }

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      }
    }
  }
  for (SpatialIndex j{0}; j < _ny; ++j) {
    for (SpatialIndex k{1}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex k{0}; k < _nz; ++k) {
    for (SpatialIndex i{1}; i < _nx; ++i) {
      for (SpatialIndex j{1}; j < _ny; ++j) {
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
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
    e->update();
  }
}

void Simulation::updateEWithDispersiveMaterial() {
  auto& ex{getEx()};
  auto& ey{getEy()};
  auto& ez{getEz()};
  auto& hx{getHx()};
  auto& hy{getHy()};
  auto& hz{getHz()};
  const auto& cexe{getCexe()};
  const auto& cexhy{getCexhy()};
  const auto& cexhz{getCexhz()};
  const auto& ceye{getCeye()};
  const auto& ceyhz{getCeyhz()};
  const auto& ceyhx{getCeyhx()};
  const auto& ceze{getCeze()};
  const auto& cezhx{getCezhx()};
  const auto& cezhy{getCezhy()};

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{1}; j < _ny; ++j) {
      for (SpatialIndex k{1}; k < _nz; ++k) {
        auto material_index{_grid_space->getGridMaterialIndex(i, j, k)};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEx(i, j, k);
          continue;
        }
        ex(i, j, k) = cexe(i, j, k) * ex(i, j, k) +
                      cexhy(i, j, k) * (hy(i, j, k) - hy(i, j, k - 1)) +
                      cexhz(i, j, k) * (hz(i, j, k) - hz(i, j - 1, k));
      }
    }
  }
  for (SpatialIndex j{0}; j < _ny; ++j) {
    for (SpatialIndex k{1}; k < _nz; ++k) {
      for (SpatialIndex i{1}; i < _nx; ++i) {
        auto material_index{_grid_space->getGridMaterialIndex(i, j, k)};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEy(i, j, k);
          continue;
        }
        ey(i, j, k) = ceye(i, j, k) * ey(i, j, k) +
                      ceyhx(i, j, k) * (hx(i, j, k) - hx(i, j, k - 1)) +
                      ceyhz(i, j, k) * (hz(i, j, k) - hz(i - 1, j, k));
      }
    }
  }
  for (SpatialIndex k{0}; k < _nz; ++k) {
    for (SpatialIndex i{1}; i < _nx; ++i) {
      for (SpatialIndex j{1}; j < _ny; ++j) {
        auto material_index{_grid_space->getGridMaterialIndex(i, j, k)};
        if (_is_exist_dispersive_material_array(material_index)) {
          _objects[material_index]->updateEz(i, j, k);
          continue;
        }
        ez(i, j, k) = ceze(i, j, k) * ez(i, j, k) +
                      cezhx(i, j, k) * (hx(i, j, k) - hx(i, j - 1, k)) +
                      cezhy(i, j, k) * (hy(i, j, k) - hy(i - 1, j, k));
      }
    }
  }
}

void Simulation::correctE() {
  for (auto&& e : _lumped_elements) {
    e->correctE();
  }

  if (_tfsf != nullptr) {
    _tfsf->updateE();
  }
}

void Simulation::correctH() {
  for (auto&& e : _lumped_elements) {
    e->correctH();
  }

  if (_tfsf != nullptr) {
    _tfsf->updateH();
  }
}

}  // namespace xfdtd
