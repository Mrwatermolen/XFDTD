#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <xtensor/xadapt.hpp>


#include "boundary/boundary.h"
#include "mesh/grid_box.h"
#include "monitor/field_monitor.h"
#include "simulation/simulation.h"
#include "simulation/yee_cell.h"
#include "util/constant.h"
#include "util/float_compare.h"
#include "util/type_define.h"

namespace xfdtd {

void Simulation::init() {
  _dt = _cfl /
        (constant::C_0 *
         std::sqrt(1.0 / (_dx * _dx) + 1.0 / (_dy * _dy) + 1.0 / (_dz * _dz)));
  _current_time_step = 0;

  for (size_t i = 0; i < _time_steps; ++i) {
    _time_array.emplace_back((i + 0.5) * _dt);
  }

  initMaterialGrid();
  initTFSF();
  initNFFFT();
  initUpdateCoefficient();
  initBoundaryCondition();
  initMonitor();
}

void Simulation::initMaterialGrid() {
  calculateDomainSize();
  gridSimulationSpace();
  allocateArray();
  _is_exist_dispersive_material_array.resize({_objects.size()});
  _is_exist_dispersive_material_array.fill(false);
  int counter_tmp{0};
  for (auto&& e : _objects) {
    e->init(_dt, _dx, getEMFInstance());
    if (e->isDispersion()) {
      _is_exist_dispersive_material = true;
      _is_exist_dispersive_material_array(counter_tmp) = true;
    }
    ++counter_tmp;
  }
  if (_is_exist_dispersive_material) {
    _emf->allocateExPrev(_nx, _ny + 1, _nz + 1);
    _emf->allocateEyPrev(_nx + 1, _ny, _nz + 1);
    _emf->allocateEzPrev(_nx + 1, _ny + 1, _nz);
    _emf->allocateJx(_nx, _ny + 1, _nz + 1);
    _emf->allocateJy(_nx + 1, _ny, _nz + 1);
    _emf->allocateJz(_nx + 1, _ny + 1, _nz);
    _emf->allocateJxPrev(_nx, _ny + 1, _nz + 1);
    _emf->allocateJyPrev(_nx + 1, _ny, _nz + 1);
    _emf->allocateJzPrev(_nx + 1, _ny + 1, _nz);
  }
}

void Simulation::initTFSF() {
  if (_tfsf == nullptr) {
    return;
  }
  auto [x, y, z]{_tfsf->getDistance()};
  _tfsf->setEMFInstance(getEMFInstance());
  _tfsf->init(_dx, _dy, _dz, _dt,
              std::make_unique<GridBox>(x, y, z, _nx - 2 * x, _ny - 2 * y,
                                        _nz - 2 * z));
}

void Simulation::initNFFFT() {
  if (_nffft == nullptr) {
    return;
  }
  auto [x, y, z]{_nffft->getDistance()};
  if (_nz == 1) {
    _nffft->init(
        std::make_unique<GridBox>(x, y, 0, _nx - 2 * x, _ny - 2 * y, 1),
        getEMFInstance(), _time_steps, _dt, _dx, _dy, 1);
    return;
  }
  _nffft->init(
      std::make_unique<GridBox>(x, y, z, _nx - 2 * x, _ny - 2 * y, _nz - 2 * z),
      getEMFInstance(), _time_steps, _dt, _dx, _dy, _dz);
}

void Simulation::initUpdateCoefficient() {
  int material_counter{0};
  for (auto&& e : _objects) {
    xt::xarray<bool> mask{xt::make_lambda_xfunction(
        [&e, &material_counter](std::shared_ptr<YeeCell>& a) {
          auto flag{e->isPointInside(a->getCenter())};
          if (flag) {
            a->setMaterialIndex(material_counter);
          }
          return flag;
        },
        _grid_space)};
    ++material_counter;
    e->correctCece(mask, _cexe, _dt);
    e->correctChch(mask, _chxh, _dt);
    e->correctCecha(mask, _cexhy, _dz, _dt);
    e->correctChcea(mask, _chxey, _dz, _dt);
    e->correctCechb(mask, _cexhz, _dy, _dt);
    e->correctChceb(mask, _chxez, _dy, _dt);

    e->correctCece(mask, _ceye, _dt);
    e->correctChch(mask, _chyh, _dt);
    e->correctCecha(mask, _ceyhz, _dx, _dt);
    e->correctChcea(mask, _chyez, _dx, _dt);
    e->correctCechb(mask, _ceyhx, _dz, _dt);
    e->correctChceb(mask, _chyex, _dz, _dt);

    e->correctCece(mask, _ceze, _dt);
    e->correctChch(mask, _chzh, _dt);
    e->correctCecha(mask, _cezhx, _dy, _dt);
    e->correctChcea(mask, _chzex, _dy, _dt);
    e->correctCechb(mask, _cezhy, _dx, _dt);
    e->correctChceb(mask, _chzey, _dx, _dt);
  }
}

void Simulation::initBoundaryCondition() {
  for (auto& e : _boundaries) {
    e->init(this);
  }
}

void Simulation::initMonitor() {
  for (auto&& e : _monitors) {
    e->setEMFInstance(getEMFInstance());
    const auto& shape{e->getShape()};
    if (shape == nullptr) {
      continue;
    }
    YeeCellArray temp;
    for (const auto& ee : _grid_space) {
      if (shape->isPointInside(ee->getCenter())) {
        temp.emplace_back(ee);
      }
    }
    e->setYeeCells(std::move(temp));
  }
}

void Simulation::calculateDomainSize() {
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::min();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::min();
  double min_z = std::numeric_limits<double>::max();
  double max_z = std::numeric_limits<double>::min();

  if (_objects.empty()) {
    throw std::runtime_error("no object in simulation");
  }

  for (const auto& e : _objects) {
    auto temp{e->getWrappedBox()};
    auto box{dynamic_cast<Cube*>(temp.get())};
    if (box == nullptr) {
      continue;
    }

    auto t{box->getXmin()};
    if (isLessOrEqual(t, min_x, constant::TOLERABLE_EPSILON)) {
      min_x = t;
    }
    t = box->getXmax();
    if (isGreaterOrEqual(t, max_x, constant::TOLERABLE_EPSILON)) {
      max_x = t;
    }
    t = box->getYmin();
    if (isLessOrEqual(t, min_y, constant::TOLERABLE_EPSILON)) {
      min_y = t;
    }
    t = box->getYmax();
    if (isGreaterOrEqual(t, max_y, constant::TOLERABLE_EPSILON)) {
      max_y = t;
    }
    t = box->getZmin();
    if (isLessOrEqual(t, min_z, constant::TOLERABLE_EPSILON)) {
      min_z = t;
    }
    t = box->getZmax();
    if (isGreaterOrEqual(t, max_z, constant::TOLERABLE_EPSILON)) {
      max_z = t;
    }
  }

  for (auto& e : _boundaries) {
    auto ori{e->getOrientation()};
    auto len{e->getSize()};
    if (ori == Orientation::XN) {
      min_x -= _dx * len;
      continue;
    }
    if (ori == Orientation::YN) {
      min_y -= _dy * len;
      continue;
    }
    if (ori == Orientation::ZN) {
      min_z -= _dz * len;
      continue;
    }
    if (ori == Orientation::XP) {
      max_x += _dx * len;
    }
    if (ori == Orientation::YP) {
      max_y += _dy * len;
    }
    if (ori == Orientation::ZP) {
      max_z += _dz * len;
    }
  }

  _nx = std::round((max_x - min_x) / _dx);
  _ny = std::round((max_y - min_y) / _dy);
  _nz = std::round((max_z - min_z) / _dz);
  double size_x = _nx * _dx;
  double size_y = _ny * _dy;
  double size_z = _nz * _dz;
  // IMPORTANT: There is a memory error while running the program with other
  // platforms.
  auto origin_point{PointVector{min_x, min_y, min_z}};
  auto box_size{PointVector{size_x, size_y, size_z}};
  _simulation_box =
      std::make_unique<Cube>(std::move(origin_point), std::move(box_size));
  if (_nx == 0) {
    _nx = 1;
    _dx = 1;
  }
  if (_ny == 0) {
    _ny = 1;
    _dy = 1;
  }
  if (_nz == 0) {
    _nz = 1;
    _dz = 1;
  }
}

void Simulation::gridSimulationSpace() {
  _grid_space.resize({static_cast<size_t>(_nx), static_cast<size_t>(_ny),
                      static_cast<size_t>(_nz)});
  auto min_x{_simulation_box->getXmin()};
  auto min_y{_simulation_box->getYmin()};
  auto min_z{_simulation_box->getZmin()};
  if (_nx == 1 && _ny == 1) {
    for (SpatialIndex k{0}; k < _nz; ++k) {
      _grid_space(0, 0, k) =
          std::make_shared<YeeCell>(-1, -1, k, _dz, min_x, min_y, min_z);
    }
  } else if (_nz == 1) {
    for (SpatialIndex i{0}; i < _nx; ++i) {
      for (SpatialIndex j{0}; j < _ny; ++j) {
        _grid_space(i, j, 0) =
            std::make_shared<YeeCell>(i, j, -1, _dx, min_x, min_y, min_z);
      }
    }
  } else {
    for (SpatialIndex i{0}; i < _nx; ++i) {
      for (SpatialIndex j{0}; j < _ny; ++j) {
        for (SpatialIndex k{0}; k < _nz; ++k) {
          _grid_space(i, j, k) =
              std::make_shared<YeeCell>(i, j, k, _dx, min_x, min_y, min_z);
        }
      }
    }
  }

}

void Simulation::allocateArray() {
  auto min_sigma{std::numeric_limits<double>::epsilon() / 1000.0};
  allocateEx(_nx, _ny + 1, _nz + 1);
  allocateEy(_nx + 1, _ny, _nz + 1);
  allocateEz(_nx + 1, _ny + 1, _nz);
  allocateHx(_nx + 1, _ny, _nz);
  allocateHy(_nx, _ny + 1, _nz);
  allocateHz(_nx, _ny, _nz + 1);

  _eps_x = allocateDoubleArray3D(_nx, _ny, _nz, constant::EPSILON_0);
  _sigma_e_x = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _mu_x = allocateDoubleArray3D(_nx, _ny, _nz, constant::MU_0);
  _sigma_m_x = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _eps_y = allocateDoubleArray3D(_nx, _ny, _nz, constant::EPSILON_0);
  _sigma_e_y = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _mu_y = allocateDoubleArray3D(_nx, _ny, _nz, constant::MU_0);
  _sigma_m_y = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _eps_z = allocateDoubleArray3D(_nx, _ny, _nz, constant::EPSILON_0);
  _sigma_e_z = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _mu_z = allocateDoubleArray3D(_nx, _ny, _nz, constant::MU_0);
  _sigma_m_z = allocateDoubleArray3D(_nx, _ny, _nz, min_sigma);

  _cexe = (2 * _eps_x - _dt * _sigma_e_x) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhz = (2 * _dt / _dy) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhy = -(2 * _dt / _dz) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexje = -(2 * _dt) / (2 * _eps_x + _dt * _sigma_e_x);

  _ceye = (2 * _eps_y - _dt * _sigma_e_y) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhx = (2 * _dt / _dz) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhz = -(2 * _dt / _dx) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyje = -(2 * _dt) / (2 * _eps_y + _dt * _sigma_e_y);

  _ceze = (2 * _eps_z - _dt * _sigma_e_z) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhy = (2 * _dt / _dx) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhx = -(2 * _dt / _dy) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezje = -(2 * _dt) / (2 * _eps_z + _dt * _sigma_e_z);

  _chxh = (2 * _mu_x - _dt * _sigma_m_x) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxey = (2 * _dt / _dz) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxez = -(2 * _dt / _dy) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxjm = -(2 * _dt) / (2 * _mu_x + _dt * _sigma_m_x);

  _chyh = (2 * _mu_y - _dt * _sigma_m_y) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyez = (2 * _dt / _dx) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyex = -(2 * _dt / _dz) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyjm = -(2 * _dt) / (2 * _mu_y + _dt * _sigma_m_y);

  _chzh = (2 * _mu_z - _dt * _sigma_m_z) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzex = (2 * _dt / _dy) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzey = -(2 * _dt / _dx) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzjm = -(2 * _dt) / (2 * _mu_z + _dt * _sigma_m_z);
}

}  // namespace xfdtd
