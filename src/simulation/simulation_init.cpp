#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <xtensor/xadapt.hpp>

#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
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
  initSource();
  initTFSF();
  initNFFFT();
  initUpdateCoefficient();
  initBondaryCondition();
  initMonitor();
}

void Simulation::initMaterialGrid() {
  caculateDomainSize();
  gridSimualtionSpace();
  allocateArray();
  caculateMaterialComponent();
}

void Simulation::initSource() {
  for (const auto& e : _sources) {
    e->init(_time_array);
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
  _nffft->init(
      std::make_unique<GridBox>(x, y, z, _nx - 2 * x, _ny - 2 * y, _nz - 2 * z),
      getEMFInstance(), _time_steps, _dt, _dx, _dy, _dz);
}

void Simulation::initUpdateCoefficient() {
  _cexe = (2 * _eps_x - _dt * _sigma_e_x) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhz = (2 * _dt / _dy) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhy = -(2 * _dt / _dz) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexje = -(2 * _dt) / (2 * _eps_x + _dt * _sigma_e_x);
  // _cexe.fill(1);
  // _cexhz.fill(_dt / (_dx * constant::EPSILON_0));
  // _cexhy.fill(-_dt / (_dx * constant::EPSILON_0));

  // _ey
  _ceye = (2 * _eps_y - _dt * _sigma_e_y) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhx = (2 * _dt / _dz) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhz = -(2 * _dt / _dx) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyje = -(2 * _dt) / (2 * _eps_y + _dt * _sigma_e_y);
  // _ceye.fill(1);
  // _ceyhx.fill(_dt / (_dz * constant::EPSILON_0));
  // _ceyhz.fill(-_dt / (_dz * constant::EPSILON_0));

  // _ez
  _ceze = (2 * _eps_z - _dt * _sigma_e_z) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhy = (2 * _dt / _dx) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhx = -(2 * _dt / _dy) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezje = -(2 * _dt) / (2 * _eps_z + _dt * _sigma_e_z);
  // _ceze.fill(1);
  // _cezhy.fill(_dt / (_dx * constant::EPSILON_0));
  // _cezhx.fill(-_dt / (_dx * constant::EPSILON_0));

  _chxh = (2 * _mu_x - _dt * _sigma_m_x) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxey = (2 * _dt / _dz) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxez = -(2 * _dt / _dy) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxjm = -(2 * _dt) / (2 * _mu_x + _dt * _sigma_m_x);
  // _chxh.fill(1);
  // _chxey.fill(_dt / (_dx * constant::MU_0));
  // _chxez.fill(-_dt / (_dx * constant::MU_0));

  // _hy
  _chyh = (2 * _mu_y - _dt * _sigma_m_y) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyez = (2 * _dt / _dx) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyex = -(2 * _dt / _dz) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyjm = -(2 * _dt) / (2 * _mu_y + _dt * _sigma_m_y);
  // _chyh.fill(1);
  // _chyez.fill(_dt / (_dx * constant::MU_0));
  // _chyex.fill(-_dt / (_dx * constant::MU_0));

  // _hz
  _chzh = (2 * _mu_z - _dt * _sigma_m_z) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzex = (2 * _dt / _dy) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzey = -(2 * _dt / _dx) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzjm = -(2 * _dt) / (2 * _mu_z + _dt * _sigma_m_z);
  // _chzh.fill(1);
  // _chzex.fill(_dt / (_dy * constant::MU_0));
  // _chzey.fill(-_dt / (_dx * constant::MU_0));
}

void Simulation::initBondaryCondition() {
  for (auto& e : _boundaries) {
    e->setEMFInstance(getEMFInstance());
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
    for (const auto& e : _grid_space) {
      if (shape->isPointInside(e->getCenter())) {
        temp.emplace_back(e);
      }
    }
    e->setYeeCells(std::move(temp));
  }
}

void Simulation::caculateDomainSize() {
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
    auto tmep{std::move(e->getWrappedBox())};
    auto box{dynamic_cast<Cube*>(tmep.get())};
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
    // tmep.release();  // For Debug
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
  max_x = min_x + size_x;
  max_y = min_y + size_y;
  max_z = min_z + size_z;
  _simulation_box = std::make_unique<Cube>(
      PointVector{min_x, min_y, min_z},
      PointVector{max_x - min_x, max_y - min_y, max_z - min_z});
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

void Simulation::gridSimualtionSpace() {
  auto total_grid{_nx * _ny * _nz};

  auto min_x{_simulation_box->getXmin()};
  auto min_y{_simulation_box->getYmin()};
  auto min_z{_simulation_box->getZmin()};
  if (_nx == 1 && _ny == 1) {
    for (SpatialIndex k{0}; k < _nz; ++k) {
      auto index{0 * _ny * _nz + 0 * _nz + k};
      _grid_space.emplace_back(
          std::make_shared<YeeCell>(PointVector{min_x, min_y, min_z + k * _dz},
                                    PointVector{0, 0, _dz}, -1, 0, 0, k));
    }
  } else if (_nz == 1) {
    for (SpatialIndex i{0}; i < _nx; ++i) {
      for (SpatialIndex j{0}; j < _ny; ++j) {
        auto index{i * _ny + j};
        _grid_space.emplace_back(std::make_shared<YeeCell>(
            PointVector{min_x + i * _dx, min_y + j * _dy, 0},
            PointVector{_dx, _dy, 0}, -1, i, j, 0));
      }
    }
  } else {
    for (SpatialIndex i{0}; i < _nx; ++i) {
      for (SpatialIndex j{0}; j < _ny; ++j) {
        for (SpatialIndex k{0}; k < _nz; ++k) {
          auto index{i * _ny * _nz + j * _nz + k};
          _grid_space.emplace_back(std::make_shared<YeeCell>(
              PointVector{min_x + i * _dx, min_y + j * _dy, min_z + k * _dz},
              PointVector{_dx, _dy, _dz}, -1, i, j, k));
        }
      }
    }
  }

  // 为每个格子设置材料
  // 允许材料覆盖
  for (auto&& c : _grid_space) {
    c->setMaterialIndex(0);
    int counter = 0;

    for (const auto& e : _objects) {
      if (!e->isPointInside(c->getCenter())) {
        ++counter;
        continue;
      }
      c->setMaterialIndex(counter);
      ++counter;
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

  _eps_x = allocateDoubleArray3D(_nx, _ny + 1, _nz + 1, constant::EPSILON_0);
  _sigma_e_x = allocateDoubleArray3D(_nx, _ny + 1, _nz + 1, min_sigma);

  _mu_x = allocateDoubleArray3D(_nx + 1, _ny, _nz, constant::MU_0);
  _sigma_m_x = allocateDoubleArray3D(_nx + 1, _ny, _nz, min_sigma);

  _eps_y = allocateDoubleArray3D(_nx + 1, _ny, _nz + 1, constant::EPSILON_0);
  _sigma_e_y = allocateDoubleArray3D(_nx + 1, _ny, _nz + 1, min_sigma);

  _mu_y = allocateDoubleArray3D(_nx, _ny + 1, _nz, constant::MU_0);
  _sigma_m_y = allocateDoubleArray3D(_nx, _ny + 1, _nz, min_sigma);

  _eps_z = allocateDoubleArray3D(_nx + 1, _ny + 1, _nz, constant::EPSILON_0);
  _sigma_e_z = allocateDoubleArray3D(_nx + 1, _ny + 1, _nz, min_sigma);

  _mu_z = allocateDoubleArray3D(_nx, _ny, _nz + 1, constant::MU_0);
  _sigma_m_z = allocateDoubleArray3D(_nx, _ny, _nz + 1, min_sigma);
}

void Simulation::caculateMaterialComponent() {
  auto min_sigma{std::numeric_limits<double>::epsilon() / 1000.0};

  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        auto material_index{
            _grid_space[i * _ny * _nz + j * _nz + k]->getMaterialIndex()};
        if (material_index == -1) {
          std::cerr << "Error: material index is -1" << std::endl;
          exit(-1);
        }
        auto [eps, mu, sigma_e, sigma_m] =
            _objects[material_index]->getElectromagneticProperties();
        // std::cout << "eps: " << eps << " mu: " << mu << " sigma_e: " <<
        // sigma_e
        //           << " sigma_m: " << sigma_m << std::endl;
        if (isLessOrEqual(sigma_e, min_sigma, constant::TOLERABLE_EPSILON)) {
          sigma_e = min_sigma;
        }
        if (isLessOrEqual(sigma_m, min_sigma, constant::TOLERABLE_EPSILON)) {
          sigma_m = min_sigma;
        }
        _eps_x(i, j, k) = eps;
        _sigma_e_x(i, j, k) = sigma_e;
        _mu_x(i, j, k) = mu;
        _sigma_m_x(i, j, k) = sigma_m;
        _eps_y(i, j, k) = eps;
        _mu_y(i, j, k) = mu;
        _sigma_e_y(i, j, k) = sigma_e;
        _sigma_m_y(i, j, k) = sigma_m;
        _eps_z(i, j, k) = eps;
        _mu_z(i, j, k) = mu;
        _sigma_e_z(i, j, k) = sigma_e;
        _sigma_m_z(i, j, k) = sigma_m;
      }
    }
  }
  // TODO(franzero): 注意这里忘记给边界赋值了
}

}  // namespace xfdtd
