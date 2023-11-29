#include <memory>
#include <xtensor/xadapt.hpp>

#include "shape/shape.h"
#include "simulation/simulation.h"
#include "util/type_define.h"

namespace xfdtd {

void Simulation::init() {
  initGridSpace();
  initFDTDBasicCoff();
  initEMInstance();
  initObject();
  initTFSF();
  initNFFFT();
  initLumpedElement();
  initUpdateCoefficient();
  initBoundaryCondition();
  initMonitor();
  initNetwork();
}

void Simulation::initFDTDBasicCoff() {
  _fdtd_basic_coff->init(_grid_space.get(), _total_time_steps);
}

void Simulation::initEMInstance() {
  _emf->allocateEx(_nx, _ny + 1, _nz + 1);
  _emf->allocateEy(_nx + 1, _ny, _nz + 1);
  _emf->allocateEz(_nx + 1, _ny + 1, _nz);
  _emf->allocateHx(_nx + 1, _ny, _nz);
  _emf->allocateHy(_nx, _ny + 1, _nz);
  _emf->allocateHz(_nx, _ny, _nz + 1);
}

void Simulation::initObject() {
  _is_exist_dispersive_material_array.resize({_objects.size()});
  _is_exist_dispersive_material_array.fill(false);
  int counter_tmp{0};
  for (auto&& e : _objects) {
    e->init(counter_tmp, _fdtd_basic_coff, _grid_space, _emf);
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
  _tfsf->init(_grid_space, _fdtd_basic_coff, _emf);
}

void Simulation::initNFFFT() {
  if (_nffft == nullptr) {
    return;
  }
  _nffft->init(_grid_space, _fdtd_basic_coff, _emf);
}

void Simulation::initLumpedElement() {
  for (auto& e : _lumped_elements) {
    e->init(_grid_space, _fdtd_basic_coff, _emf);
  }
}

void Simulation::initUpdateCoefficient() {
  _fdtd_basic_coff->initCoff(_grid_space.get());
  for (auto&& e : _objects) {
    e->correctFDTDCoff();
  }

  for (auto&& e : _lumped_elements) {
    e->correctFDTDCoff();
  }
}

void Simulation::initBoundaryCondition() {
  for (auto& e : _boundaries) {
    e->init(_emf, _fdtd_basic_coff, _grid_space, nullptr);
  }
}

void Simulation::initMonitor() {
  for (auto&& e : _monitors) {
    e->init(_fdtd_basic_coff, _grid_space, _emf);
  }
}

void Simulation::initGridSpace() {
  std::vector<std::unique_ptr<Shape>> shapes;
  for (const auto& object : _objects) {
    _grid_space->calculateSpaceSize(object->getWrappedBox().get());
  }

  auto base_dx{_grid_space->getGridBaseSizeX()};
  auto base_dy{_grid_space->getGridBaseSizeY()};
  auto base_dz{_grid_space->getGridBaseSizeZ()};
  for (const auto& boundary : _boundaries) {
    auto ori{boundary->getOrientation()};
    if (ori == Orientation::XN) {
      _grid_space->extendBoundaryXN(boundary->getSize() * base_dx);
      continue;
    }
    if (ori == Orientation::XP) {
      _grid_space->extendBoundaryXP(boundary->getSize() * base_dx);
      continue;
    }
    if (ori == Orientation::YN) {
      _grid_space->extendBoundaryYN(boundary->getSize() * base_dy);
      continue;
    }
    if (ori == Orientation::YP) {
      _grid_space->extendBoundaryYP(boundary->getSize() * base_dy);
      continue;
    }
    if (ori == Orientation::ZN) {
      _grid_space->extendBoundaryZN(boundary->getSize() * base_dz);
      continue;
    }
    if (ori == Orientation::ZP) {
      _grid_space->extendBoundaryZP(boundary->getSize() * base_dz);
      continue;
    }
  }
  _grid_space->generateUniformGridSpace();
  _nx = _grid_space->getGridNumX();
  _ny = _grid_space->getGridNumY();
  _nz = _grid_space->getGridNumZ();
  _grid_space->tellMeOk();
}

void Simulation::initNetwork() {
  if (_network == nullptr) {
    return;
  }
  _network->init(_fdtd_basic_coff, _grid_space, _emf);
}

}  // namespace xfdtd
