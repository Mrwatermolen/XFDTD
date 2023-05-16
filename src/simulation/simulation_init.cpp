#include <cmath>
#include <exception>
#include <iostream>
#include <memory>

#include "boundary/boundary.h"
#include "boundary/perfect_match_layer.h"
#include "util/constant.h"
#include "simulation/simulation.h"
#include "simulation/yee_cell.h"
#include "util/float_compare.h"

#ifdef DEBUG

#include <iostream>

#endif  // DEBUG

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
  initUpdateCoefficient();
  initBondaryCondition();
  initOutputParameter();

  _ex.setZero();
  _ey.setZero();
  _ez.setZero();
  _hx.setZero();
  _hy.setZero();
  _hz.setZero();
}

void Simulation::initMaterialGrid() {
  caculateDomainSize();
  gridSimualtionSpace();
  caculateMaterialComponent();
}

void Simulation::initSource() {
  for (const auto& e : _sources) {
    e->init(_time_array);
  }
}

void Simulation::initUpdateCoefficient() {
  _cexe.resize(_nx, _ny + 1, _nz + 1);
  _cexhy.resize(_nx, _ny + 1, _nz + 1);
  _cexhz.resize(_nx, _ny + 1, _nz + 1);
  _ceye.resize(_nx + 1, _ny, _nz + 1);
  _ceyhx.resize(_nx + 1, _ny, _nz + 1);
  _ceyhz.resize(_nx + 1, _ny, _nz + 1);
  _ceze.resize(_nx + 1, _ny + 1, _nz);
  _cezhx.resize(_nx + 1, _ny + 1, _nz);
  _cezhy.resize(_nx + 1, _ny + 1, _nz);

  _chxh.resize(_nx, _ny + 1, _nz);
  _chxey.resize(_nx, _ny + 1, _nz);
  _chxez.resize(_nx, _ny + 1, _nz);
  _chyh.resize(_nx + 1, _ny, _nz);
  _chyez.resize(_nx + 1, _ny, _nz);
  _chyex.resize(_nx + 1, _ny, _nz);
  _chzh.resize(_nx, _ny, _nz);
  _chzex.resize(_nx, _ny, _nz + 1);
  _chzey.resize(_nx, _ny, _nz + 1);

  _cexe = (2 * _eps_x - _dt * _sigma_e_x) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhz = (2 * _dt / _dy) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexhy = -(2 * _dt / _dz) / (2 * _eps_x + _dt * _sigma_e_x);
  _cexje = -(2 * _dt) / (2 * _eps_x + _dt * _sigma_e_x);

  // _ey
  _ceye = (2 * _eps_y - _dt * _sigma_e_y) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhx = (2 * _dt / _dz) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyhz = -(2 * _dt / _dx) / (2 * _eps_y + _dt * _sigma_e_y);
  _ceyje = -(2 * _dt) / (2 * _eps_y + _dt * _sigma_e_y);

  // _ez
  _ceze = (2 * _eps_z - _dt * _sigma_e_z) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhx = -(2 * _dt / _dy) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezhy = (2 * _dt / _dx) / (2 * _eps_z + _dt * _sigma_e_z);
  _cezje = -(2 * _dt) / (2 * _eps_z + _dt * _sigma_e_z);

  _chxh = (2 * _mu_x - _dt * _sigma_m_x) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxey = (2 * _dt / _dz) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxez = -(2 * _dt / _dy) / (2 * _mu_x + _dt * _sigma_m_x);
  _chxjm = -(2 * _dt) / (2 * _mu_x + _dt * _sigma_m_x);

  // _hy
  _chyh = (2 * _mu_y - _dt * _sigma_m_y) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyez = (2 * _dt / _dx) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyex = -(2 * _dt / _dz) / (2 * _mu_y + _dt * _sigma_m_y);
  _chyjm = -(2 * _dt) / (2 * _mu_y + _dt * _sigma_m_y);

  // _hz
  _chzh = (2 * _mu_z - _dt * _sigma_m_z) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzex = (2 * _dt / _dy) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzey = -(2 * _dt / _dx) / (2 * _mu_z + _dt * _sigma_m_z);
  _chzjm = -(2 * _dt) / (2 * _mu_z + _dt * _sigma_m_z);
}

void Simulation::initBondaryCondition() {
  for (auto& e : _boundaries) {
    auto cpml{std::dynamic_pointer_cast<PML>(e)};
    if (cpml == nullptr) {
      throw std::exception();
    }
    auto ori{e->getOrientation()};
    if (ori == Orientation::XN) {
      cpml->init(_dx, _dt, 0, _ny, _nz, _ceyhz, _cezhy, _chyez, _chzey);
      continue;
    }
    if (ori == Orientation::YN) {
      cpml->init(_dy, _dt, 0, _nz, _nx, _cezhx, _cexhz, _chzex, _chxez);
    }
    if (ori == Orientation::ZN) {
      cpml->init(_dz, _dt, 0, _nx, _ny, _cexhy, _ceyhx, _chxey, _chyex);
    }
    if (ori == Orientation::XP) {
      cpml->init(_dx, _dt, _nx - e->getSize(), _ny, _nz, _ceyhz, _cezhy, _chyez,
                 _chzey);
    }
    if (ori == Orientation::YP) {
      cpml->init(_dy, _dt, _ny - e->getSize(), _nz, _nx, _cezhx, _cexhz, _chzex,
                 _chxez);
    }
    if (ori == Orientation::ZP) {
      cpml->init(_dz, _dt, _nz - e->getSize(), _nx, _ny, _cexhy, _ceyhx, _chxey,
                 _chyex);
    }
  }
}

void Simulation::initOutputParameter() {}

void Simulation::caculateDomainSize() {
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::min();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::min();
  double min_z = std::numeric_limits<double>::max();
  double max_z = std::numeric_limits<double>::min();

  for (const auto& e : _objects) {
    std::shared_ptr<Shape> tmep{std::move(e->getWrappedBox())};
    auto box{std::dynamic_pointer_cast<Cube>(tmep)};
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

  auto nx = std::round((max_x - min_x) / _dx);
  auto ny = std::round((max_y - min_y) / _dy);
  auto nz = std::round((max_z - min_z) / _dz);
  double size_x = nx * _dx;
  double size_y = ny * _dy;
  double size_z = nz * _dz;

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
  _simulation_box = std::make_unique<Cube>(
      Eigen::Vector3d{min_x, min_y, min_z},
      Eigen::Vector3d{max_x - min_x, max_y - min_y, max_z - min_z});
}

void Simulation::gridSimualtionSpace() {
  _nx = std::round(_simulation_box->getSize().x() / _dx);
  _ny = std::round(_simulation_box->getSize().y() / _dy);
  _nz = std::round(_simulation_box->getSize().z() / _dz);
  if (_nx == 0) {
    _nx = 1;
  }
  if (_ny == 0) {
    _ny = 1;
  }
  if (_nz == 0) {
    _nz = 1;
  }
  // check boundary

  auto total_grid{_nx * _ny * _nz};

  auto min_sigma{std::numeric_limits<double>::epsilon() / 1000.0};
  _ex = EFTA(_nx, _ny + 1, _nz + 1);
  _eps_x = EFTA(_nx, _ny + 1, _nz + 1);
  _eps_x.setConstant(constant::EPS_0);
  _sigma_e_x = EFTA(_nx, _ny + 1, _nz + 1);
  _sigma_e_x.setConstant(min_sigma);

  _hx = EFTA(_nx + 1, _ny, _nz);
  _mu_x = EFTA(_nx + 1, _ny, _nz);
  _mu_x.setConstant(constant::MU_0);
  _sigma_m_x = EFTA(_nx + 1, _ny, _nz);
  _sigma_m_x = _sigma_m_x.setConstant(min_sigma);

  _ey = EFTA(_nx + 1, _ny, _nz + 1);
  _eps_y = EFTA(_nx + 1, _ny, _nz + 1);
  _eps_y.setConstant(constant::EPS_0);
  _sigma_e_y = EFTA(_nx + 1, _ny, _nz + 1);
  _sigma_e_y.setConstant(min_sigma);

  _hy = EFTA(_nx, _ny + 1, _nz);
  _mu_y = EFTA(_nx, _ny + 1, _nz);
  _mu_y.setConstant(constant::MU_0);
  _sigma_m_y = EFTA(_nx, _ny + 1, _nz);
  _sigma_m_y.setConstant(min_sigma);

  _ez = EFTA(_nx + 1, _ny + 1, _nz);
  _eps_z = EFTA(_nx + 1, _ny + 1, _nz);
  _eps_z.setConstant(constant::EPS_0);
  _sigma_e_z = EFTA(_nx + 1, _ny + 1, _nz);
  _sigma_e_z.setConstant(min_sigma);

  _hz = EFTA(_nx, _ny, _nz + 1);
  _mu_z = EFTA(_nx, _ny, _nz + 1);
  _mu_z.setConstant(constant::MU_0);
  _sigma_m_z = EFTA(_nx, _ny, _nz + 1);
  _sigma_m_z.setConstant(min_sigma);

  auto min_x{_simulation_box->getXmin()};
  auto min_y{_simulation_box->getYmin()};
  auto min_z{_simulation_box->getZmin()};
  // 1d
  if (_nx == 1 && _ny == 1) {
    for (SpatialIndex k{0}; k < _nz; ++k) {
      auto index{0 * _ny * _nz + 0 * _nz + k};
      _grid_space.emplace_back(std::make_shared<YeeCell>(
          Eigen::Vector3d{min_x, min_y, min_z + k * _dz},
          Eigen::Vector3d{0, 0, _dz}, 0));
    }
  } else {
    for (SpatialIndex i{0}; i < _nx; ++i) {
      for (SpatialIndex j{0}; j < _ny; ++j) {
        for (SpatialIndex k{0}; k < _nz; ++k) {
          auto index{i * _ny * _nz + j * _nz + k};
          _grid_space.emplace_back(std::make_shared<YeeCell>(
              Eigen::Vector3d{min_x + i * _dx, min_y + j * _dy,
                              min_z + k * _dz},
              Eigen::Vector3d{_dx, _dy, _dz}, -1));
        }
      }
    }
  }

  // 为每个格子设置材料
  // 允许材料覆盖
  for (auto&& c : _grid_space) {
    // c->setMaterialIndex(0);
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

void Simulation::caculateMaterialComponent() {
  auto min_sigma{std::numeric_limits<double>::epsilon() / 1000.0};
  for (SpatialIndex i{0}; i < _nx; ++i) {
    for (SpatialIndex j{0}; j < _ny; ++j) {
      for (SpatialIndex k{0}; k < _nz; ++k) {
        auto [eps, mu, sigma_e, sigma_m] =
            _objects[_grid_space[i * _ny * _nz + j * _nz + k]
                         ->getMaterialIndex()]
                ->getElectromagneticProperties();
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
