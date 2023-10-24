#include "fdtd_basic_coff/fdtd_basic_coff.h"

#include <cstddef>

#include "util/constant.h"

namespace xfdtd {

FDTDBasicCoff::FDTDBasicCoff(double cfl) : _cfl{cfl} {}

void FDTDBasicCoff::init(const GridSpace* grid_space, size_t total_time_step,
                         size_t current_time_step) {
  if (grid_space == nullptr) {
    throw std::runtime_error("Grid space is not initialized.");
  }

  setTotalTimeStep(total_time_step);
  setCurrentTimeStep(current_time_step);
  if (grid_space->getGridNumZ() <= 1) {
    setDt(calculateDeltaT(_cfl, grid_space->getGridSizeMinX(),
                          grid_space->getGridSizeMinY()));
  } else {
    setDt(calculateDeltaT(_cfl, grid_space->getGridSizeMinX(),
                          grid_space->getGridSizeMinY(),
                          grid_space->getGridSizeMinZ()));
  }
  allocateArray(grid_space->getGridNumX(), grid_space->getGridNumY(),
                grid_space->getGridNumZ());
}

void FDTDBasicCoff::initCoff(const GridSpace* grid_space) {
  if (grid_space == nullptr) {
    throw std::runtime_error("Grid space is not initialized.");
  }
  // auto dx{grid_space->getGridBaseSizeX()};
  // auto dy{grid_space->getGridBaseSizeY()};
  // auto dz{grid_space->getGridBaseSizeZ()};
  // _dx = dx;
  // _dy = dy;
  // _dz = dz;
  // _cexe = (2 * _eps_x - _dt * _sigma_x) / (2 * _eps_x + _dt * _sigma_x);
  // _cexhy = -(2 * _dt / _dz) / (2 * _eps_x + _sigma_x * _dt);
  // _cexhz = (2 * _dt / _dy) / (2 * _eps_x + _sigma_x * _dt);

  // _ceye = (2 * _eps_y - _dt * _sigma_y) / (2 * _eps_y + _dt * _sigma_y);
  // _ceyhz = -(2 * _dt / _dx) / (2 * _eps_y + _sigma_y * _dt);
  // _ceyhx = (2 * _dt / _dz) / (2 * _eps_y + _sigma_y * _dt);

  // _ceze = (2 * _eps_z - _dt * _sigma_z) / (2 * _eps_z + _dt * _sigma_z);
  // _cezhx = -(2 * _dt / _dy) / (2 * _eps_z + _sigma_z * _dt);
  // _cezhy = (2 * _dt / _dx) / (2 * _eps_z + _sigma_z * _dt);

  // _chxh = (2 * _mu_x - _dt * _sigma_m_x) / (2 * _mu_x + _dt * _sigma_m_x);
  // _chxey = (2 * _dt / _dy) / (2 * _mu_x + _sigma_m_x * _dt);
  // _chxez = -(2 * _dt / _dz) / (2 * _mu_x + _sigma_m_x * _dt);

  // _chyh = (2 * _mu_y - _dt * _sigma_m_y) / (2 * _mu_y + _dt * _sigma_m_y);
  // _chyez = (2 * _dt / _dz) / (2 * _mu_y + _sigma_m_y * _dt);
  // _chyex = -(2 * _dt / _dx) / (2 * _mu_y + _sigma_m_y * _dt);

  // _chzh = (2 * _mu_z - _dt * _sigma_m_z) / (2 * _mu_z + _dt * _sigma_m_z);
  // _chzex = (2 * _dt / _dx) / (2 * _mu_z + _sigma_m_z * _dt);
  // _chzey = -(2 * _dt / _dy) / (2 * _mu_z + _sigma_m_z * _dt);
  // return;

  const auto& e_node_size_x{grid_space->getGridSizeArrayEX()};
  const auto& e_node_size_y{grid_space->getGridSizeArrayEY()};
  const auto& e_node_size_z{grid_space->getGridSizeArrayEZ()};
  const auto& h_node_size_x{grid_space->getGridSizeArrayHX()};
  const auto& h_node_size_y{grid_space->getGridSizeArrayHY()};
  const auto& h_node_size_z{grid_space->getGridSizeArrayHZ()};

  auto get_xyz{[](const auto& size_x, const auto& size_y, const auto& size_z) {
    return xt::meshgrid(size_x, size_y, size_z);
  }};

  auto calculate_coff_e{[](const auto& da, const auto& db, const auto& dc,
                           const auto& dt, const auto& eps, const auto& sigma,
                           auto&& cece, auto&& cecha, auto&& cechb) {
    cece = (2 * eps - dt * sigma) / (2 * eps + dt * sigma);
    cecha = -(2 * dt / db) / (2 * eps + sigma * dt);
    cechb = (2 * dt / da) / (2 * eps + sigma * dt);
  }};

  auto calculate_coff_h{[](const auto& da, const auto& db, const auto& dc,
                           const auto& dt, const auto& mu, const auto& sigma_m,
                           auto&& chch, auto&& chcea, auto&& chceb) {
    chch = (2 * mu - dt * sigma_m) / (2 * mu + dt * sigma_m);
    chcea = (2 * dt / db) / (2 * mu + sigma_m * dt);
    chceb = -(2 * dt / da) / (2 * mu + sigma_m * dt);
  }};

  auto e_x_size{get_xyz(e_node_size_x, h_node_size_y, h_node_size_z)};
  calculate_coff_e(std::get<1>(e_x_size), std::get<2>(e_x_size),
                   std::get<0>(e_x_size), _dt, _eps_x, _sigma_x, _cexe, _cexhy,
                   _cexhz);

  auto e_y_size{get_xyz(h_node_size_x, e_node_size_y, h_node_size_z)};
  calculate_coff_e(std::get<2>(e_y_size), std::get<0>(e_y_size),
                   std::get<1>(e_y_size), _dt, _eps_y, _sigma_y, _ceye, _ceyhz,
                   _ceyhx);

  auto e_z_size{get_xyz(h_node_size_x, h_node_size_y, e_node_size_z)};
  calculate_coff_e(std::get<0>(e_z_size), std::get<1>(e_z_size),
                   std::get<2>(e_z_size), _dt, _eps_z, _sigma_z, _ceze, _cezhx,
                   _cezhy);

  auto h_x_size{get_xyz(h_node_size_x, e_node_size_y, e_node_size_z)};
  calculate_coff_h(std::get<1>(h_x_size), std::get<2>(h_x_size),
                   std::get<0>(h_x_size), _dt, _mu_x, _sigma_m_x, _chxh, _chxey,
                   _chxez);

  auto h_y_size{get_xyz(e_node_size_x, h_node_size_y, e_node_size_z)};
  calculate_coff_h(std::get<2>(h_y_size), std::get<0>(h_y_size),
                   std::get<1>(h_y_size), _dt, _mu_y, _sigma_m_y, _chyh, _chyez,
                   _chyex);

  auto h_z_size{get_xyz(e_node_size_x, e_node_size_y, h_node_size_z)};
  calculate_coff_h(std::get<0>(h_z_size), std::get<1>(h_z_size),
                   std::get<2>(h_z_size), _dt, _mu_z, _sigma_m_z, _chzh, _chzex,
                   _chzey);
}

double FDTDBasicCoff::getCFL() const { return _cfl; }

double FDTDBasicCoff::getDt() const { return _dt; }

size_t FDTDBasicCoff::getTotalTimeStep() const { return _total_time_step; }

size_t FDTDBasicCoff::getCurrentTimeStep() const { return _current_time_step; }

const xt::xarray<double>& FDTDBasicCoff::getEpsX() const { return _eps_x; }

const xt::xarray<double>& FDTDBasicCoff::getEpsY() const { return _eps_y; }

const xt::xarray<double>& FDTDBasicCoff::getEpsZ() const { return _eps_z; }

const xt::xarray<double>& FDTDBasicCoff::getMuX() const { return _mu_x; }

const xt::xarray<double>& FDTDBasicCoff::getMuY() const { return _mu_y; }

const xt::xarray<double>& FDTDBasicCoff::getMuZ() const { return _mu_z; }

const xt::xarray<double>& FDTDBasicCoff::getSigmaX() const { return _sigma_x; }

const xt::xarray<double>& FDTDBasicCoff::getSigmaY() const { return _sigma_y; }

const xt::xarray<double>& FDTDBasicCoff::getSigmaZ() const { return _sigma_z; }

const xt::xarray<double>& FDTDBasicCoff::getSigmaMx() const {
  return _sigma_m_x;
}

const xt::xarray<double>& FDTDBasicCoff::getSigmaMy() const {
  return _sigma_m_y;
}

const xt::xarray<double>& FDTDBasicCoff::getSigmaMz() const {
  return _sigma_m_z;
}

const xt::xarray<double>& FDTDBasicCoff::getCexe() const { return _cexe; }

const xt::xarray<double>& FDTDBasicCoff::getCexhy() const { return _cexhy; }

const xt::xarray<double>& FDTDBasicCoff::getCexhz() const { return _cexhz; }

const xt::xarray<double>& FDTDBasicCoff::getCeye() const { return _ceye; }

const xt::xarray<double>& FDTDBasicCoff::getCeyhz() const { return _ceyhz; }

const xt::xarray<double>& FDTDBasicCoff::getCeyhx() const { return _ceyhx; }

const xt::xarray<double>& FDTDBasicCoff::getCeze() const { return _ceze; }

const xt::xarray<double>& FDTDBasicCoff::getCezhx() const { return _cezhx; }

const xt::xarray<double>& FDTDBasicCoff::getCezhy() const { return _cezhy; }

const xt::xarray<double>& FDTDBasicCoff::getChxh() const { return _chxh; }

const xt::xarray<double>& FDTDBasicCoff::getChxey() const { return _chxey; }

const xt::xarray<double>& FDTDBasicCoff::getChxez() const { return _chxez; }

const xt::xarray<double>& FDTDBasicCoff::getChyh() const { return _chyh; }

const xt::xarray<double>& FDTDBasicCoff::getChyez() const { return _chyez; }

const xt::xarray<double>& FDTDBasicCoff::getChyex() const { return _chyex; }

const xt::xarray<double>& FDTDBasicCoff::getChzh() const { return _chzh; }

const xt::xarray<double>& FDTDBasicCoff::getChzex() const { return _chzex; }

const xt::xarray<double>& FDTDBasicCoff::getChzey() const { return _chzey; }

xt::xarray<double>& FDTDBasicCoff::getEpsX() { return _eps_x; }

xt::xarray<double>& FDTDBasicCoff::getEpsY() { return _eps_y; }

xt::xarray<double>& FDTDBasicCoff::getEpsZ() { return _eps_z; }

xt::xarray<double>& FDTDBasicCoff::getMuX() { return _mu_x; }

xt::xarray<double>& FDTDBasicCoff::getMuY() { return _mu_y; }

xt::xarray<double>& FDTDBasicCoff::getMuZ() { return _mu_z; }

xt::xarray<double>& FDTDBasicCoff::getSigmaX() { return _sigma_x; }

xt::xarray<double>& FDTDBasicCoff::getSigmaY() { return _sigma_y; }

xt::xarray<double>& FDTDBasicCoff::getSigmaZ() { return _sigma_z; }

xt::xarray<double>& FDTDBasicCoff::getSigmaMx() { return _sigma_m_x; }

xt::xarray<double>& FDTDBasicCoff::getSigmaMy() { return _sigma_m_y; }

xt::xarray<double>& FDTDBasicCoff::getSigmaMz() { return _sigma_m_z; }

xt::xarray<double>& FDTDBasicCoff::getCexe() { return _cexe; }

xt::xarray<double>& FDTDBasicCoff::getCexhy() { return _cexhy; }

xt::xarray<double>& FDTDBasicCoff::getCexhz() { return _cexhz; }

xt::xarray<double>& FDTDBasicCoff::getCeye() { return _ceye; }

xt::xarray<double>& FDTDBasicCoff::getCeyhz() { return _ceyhz; }

xt::xarray<double>& FDTDBasicCoff::getCeyhx() { return _ceyhx; }

xt::xarray<double>& FDTDBasicCoff::getCeze() { return _ceze; }

xt::xarray<double>& FDTDBasicCoff::getCezhx() { return _cezhx; }

xt::xarray<double>& FDTDBasicCoff::getCezhy() { return _cezhy; }

xt::xarray<double>& FDTDBasicCoff::getChxh() { return _chxh; }

xt::xarray<double>& FDTDBasicCoff::getChxey() { return _chxey; }

xt::xarray<double>& FDTDBasicCoff::getChxez() { return _chxez; }

xt::xarray<double>& FDTDBasicCoff::getChyh() { return _chyh; }

xt::xarray<double>& FDTDBasicCoff::getChyez() { return _chyez; }

xt::xarray<double>& FDTDBasicCoff::getChyex() { return _chyex; }

xt::xarray<double>& FDTDBasicCoff::getChzh() { return _chzh; }

xt::xarray<double>& FDTDBasicCoff::getChzex() { return _chzex; }

xt::xarray<double>& FDTDBasicCoff::getChzey() { return _chzey; }

void FDTDBasicCoff::setCFL(double cfl) { _cfl = cfl; }

void FDTDBasicCoff::setDt(double dt) { _dt = dt; }

void FDTDBasicCoff::setTotalTimeStep(size_t total_time_step) {
  _total_time_step = total_time_step;
}

void FDTDBasicCoff::setCurrentTimeStep(size_t current_time_step) {
  _current_time_step = current_time_step;
}

void FDTDBasicCoff::allocateArray(size_t nx, size_t ny, size_t nz) {
  // _cexe.resize({nx, ny + 1, nz + 1});
  // _cexe.fill(1);
  // _cexhz.resize({nx, ny + 1, nz + 1});
  // _cexhz.fill(_dt / (_dx * constant::EPSILON_0));
  // _cexhy.resize({nx, ny + 1, nz + 1});
  // _cexhy.fill(-_dt / (_dx * constant::EPSILON_0));

  // _ceye.resize({nx + 1, ny, nz + 1});
  // _ceye.fill(1);
  // _ceyhx.resize({nx + 1, ny, nz + 1});
  // _ceyhx.fill(_dt / (_dz * constant::EPSILON_0));
  // _ceyhz.resize({nx + 1, ny, nz + 1});
  // _ceyhz.fill(-_dt / (_dz * constant::EPSILON_0));

  // _ceze.resize({nx + 1, ny + 1, nz});
  // _cezhx.resize({nx + 1, ny + 1, nz});
  // _cezhy.resize({nx + 1, ny + 1, nz});
  // _ceze.fill(1);
  // _cezhy.fill(_dt / (_dx * constant::EPSILON_0));
  // _cezhx.fill(-_dt / (_dx * constant::EPSILON_0));

  // _chxh.resize({nx + 1, ny, nz});
  // _chxey.resize({nx + 1, ny, nz});
  // _chxez.resize({nx + 1, ny, nz});
  // _chxh.fill(1);
  // _chxey.fill(_dt / (_dx * constant::MU_0));
  // _chxez.fill(-_dt / (_dx * constant::MU_0));

  // _chyh.resize({nx, ny + 1, nz});
  // _chyex.resize({nx, ny + 1, nz});
  // _chyez.resize({nx, ny + 1, nz});
  // _chyh.fill(1);
  // _chyez.fill(_dt / (_dx * constant::MU_0));
  // _chyex.fill(-_dt / (_dx * constant::MU_0));

  // _chzh.resize({nx, ny, nz + 1});
  // _chzex.resize({nx, ny, nz + 1});
  // _chzey.resize({nx, ny, nz + 1});
  // _chzh.fill(1);
  // _chzex.fill(_dt / (_dy * constant::MU_0));
  // _chzey.fill(-_dt / (_dx * constant::MU_0));

  _eps_x.resize({nx, ny + 1, nz + 1});
  _eps_x.fill(constant::EPSILON_0);
  _eps_y.resize({nx + 1, ny, nz + 1});
  _eps_y.fill(constant::EPSILON_0);
  _eps_z.resize({nx + 1, ny + 1, nz});
  _eps_z.fill(constant::EPSILON_0);

  _mu_x.resize({nx + 1, ny, nz});
  _mu_x.fill(constant::MU_0);
  _mu_y.resize({nx, ny + 1, nz});
  _mu_y.fill(constant::MU_0);
  _mu_z.resize({nx, ny, nz + 1});
  _mu_z.fill(constant::MU_0);

  _sigma_x.resize({nx, ny + 1, nz + 1});
  _sigma_x.fill(0);
  _sigma_y.resize({nx + 1, ny, nz + 1});
  _sigma_y.fill(0);
  _sigma_z.resize({nx + 1, ny + 1, nz});
  _sigma_z.fill(0);

  _sigma_m_x.resize({nx + 1, ny, nz});
  _sigma_m_x.fill(0);
  _sigma_m_y.resize({nx, ny + 1, nz});
  _sigma_m_y.fill(0);
  _sigma_m_z.resize({nx, ny, nz + 1});
  _sigma_m_z.fill(0);
}

double FDTDBasicCoff::calculateDeltaT(double cfl, double dx, double dy,
                                      double dz) {
  return cfl / (constant::C_0 *
                std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
}

double FDTDBasicCoff::calculateDeltaT(double cfl, double dx, double dy) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy)));
}

double FDTDBasicCoff::calculateDeltaT(double cfl, double dx) {
  return cfl / (constant::C_0 * std::sqrt(1.0 / (dx * dx)));
}
}  // namespace xfdtd
